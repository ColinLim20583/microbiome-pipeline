import os
import gzip
import subprocess
import itertools
import re
from functools import lru_cache
from Bio.Data.IUPACData import ambiguous_dna_values

# Define primer sets
PRIMERS_16S = [
    # 338F / 806R – V3–V4 (~460bp)
    ("ACTCCTACGGGAGGCAGCAG", "GGACTACHVGGGTWTCTAAT"),
    # 515fmod / 806rmod – V4 (~300bp)
    ("GTGYCAGCMGCCGCGGTAA", "GGACTACNVGGGTWTCTAAT"),
    # 515F / 926R – V4–V5 (~400bp)
    ("GTGYCAGCMGCCGCGGTAA", "CCGYCAATTYMTTTRAGTTT"),
    # 515F / 907R – V4–V5 (~400bp)
    ("GTGCCAGCMGCCGCGGTAA", "CCGTCAATTCMTTTRAGTTT"),
    # 341F / 785R – V3–V4 (~440bp)
    ("CCTACGGGNGGCWGCAG", "GACTACHVGGGTATCTAATCC"),
    # 799F / 1193R – V5–V7 (~400bp, reduces plant chloroplast)
    ("AACMGGATTAGATACCCKG", "ACGTCATCCCCACCTTCC"),
]

PRIMERS_ITS = [
    # ITS2-F / ITS2-R (~340bp)
    ("GTGAATCATCGARTCTTTG", "TCCTCCGCTTATTGATATGC"),
    # ITS1F / ITS2 (~310bp)
    ("CTTGGTCATTTAGAGGAAGTAA", "GCTGCGTTCTTCATCGATGC"),
    # ITS3 / ITS4 (~300bp, targets ITS2 region)
    ("GCATCGATGAAGAACGCAGC", "TCCTCCGCTTATTGATATGC"),
    # Variant of ITS2 (rare, but included)
    ("TGCGTTCTTCATCGATGC", "TCCTCCGCTTATTGATATGC"),
]

PRIMERS_AMF = [
    ("ATCAACTTTCGATGGTAGGATAGA", "GAACCCAAACACTTTGGTTTCC")   # AML1 / AML2
]

PRIMERS_18S = [
    ("CGWTAACGAACGAGACCT", "AICCATTCAATCGGTAIT"),               # FF390 / FR1
    ("GGCAAGTCTGGTGCCAG", "GACTACGACGGTATCTRATCRTCTTCG"),       # MEG18SV4F / MEG18SV4R
    ("CCAGCASCYGCGGTAATTCC", "ACTTTCGTTCTTGATYRA")              # TAReuk454FWD1 / TAReukREV3
]

@lru_cache(maxsize=None)
def expand_degenerate(primer_seq):
    # Add support for inosine (I). Biopython's ambiguous_dna_values doesn't include it.
    # Common practice: treat I as matching any base (ACGT).
    local_map = dict(ambiguous_dna_values)
    local_map["I"] = "ACGT"
    local_map["i"] = "ACGT"

    bases = [local_map[base] for base in primer_seq]
    return [''.join(p) for p in itertools.product(*bases)]

def read_first_seqs(file_path, n=1000):
    sequences = []
    open_func = gzip.open if file_path.endswith(".gz") else open
    try:
        with open_func(file_path, "rt") as f:
            for i, line in enumerate(f):
                if i % 4 == 1:
                    seq = line.strip()
                    if seq:
                        sequences.append(seq)
                    if len(sequences) >= n:
                        break
    except Exception as e:
        print(f"❌ Error reading {file_path}: {e}")
    return sequences

def run_cutadapt(fwd_file, rev_file, fwd_primers, rev_primers, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    sample_id = "_".join(os.path.basename(fwd_file).split("_")[:-1])
    out_f = os.path.join(output_dir, f"{sample_id}_R1_trimmed.fastq")
    out_r = os.path.join(output_dir, f"{sample_id}_R2_trimmed.fastq")

    cmd = ["cutadapt"]
    for p in fwd_primers:
        cmd += ["-g", p]
    for p in rev_primers:
        cmd += ["-G", p]
    cmd += ["-o", out_f, "-p", out_r, fwd_file, rev_file, "--discard-untrimmed"]

    print(f"✂️ Running Cutadapt: {sample_id}")
    try:
        subprocess.run(cmd, check=True)
        print(f"✅ Cutadapt completed for {sample_id}")
    except subprocess.CalledProcessError as e:
        print(f"❌ Cutadapt failed for {sample_id}: {e}")

def detect_primer_pair(seqs, primer_pairs):
    max_hits = 0
    best_pair = None
    for fwd, rev in primer_pairs:
        fwd_hits = sum(any(p[:12] in seq[:100] for p in expand_degenerate(fwd)) for seq in seqs)
        rev_hits = sum(any(p[:12] in seq[:100] for p in expand_degenerate(rev)) for seq in seqs)
        total = fwd_hits + rev_hits
        if total > max_hits:
            max_hits = total
            best_pair = (fwd, rev)
    return best_pair

def process_folder(folder, mode, output_base="data"):
    print(f"\n📂 Processing folder: {folder} [{mode}]")
    files = sorted([f for f in os.listdir(folder) if f.endswith(".fastq") or f.endswith(".fastq.gz")])
    paired_files = {}

    for f in files:
        match = re.match(r"(.*)_[12]", f)
        if match:
            sample = match.group(1)
            paired_files.setdefault(sample, []).append(f)

    for sample, reads in paired_files.items():
        if len(reads) != 2:
            print(f"⚠️ Skipping {sample}: expected 2 reads, found {len(reads)}")
            continue

        fwd = os.path.join(folder, sorted(reads)[0])
        rev = os.path.join(folder, sorted(reads)[1])
        seqs = read_first_seqs(fwd)

        if mode == "16S":
            best = detect_primer_pair(seqs, PRIMERS_16S)
        elif mode == "ITS":
            best = detect_primer_pair(seqs, PRIMERS_ITS)
        elif mode == "AMF":
            print(f"ℹ️ Forcing AML1/AML2 primers for {sample} (primer detection skipped).")
            best = PRIMERS_AMF[0]
        elif mode == "18S":
            best = detect_primer_pair(seqs, PRIMERS_18S)
        else:
            print(f"⚠️ Unknown mode: {mode}")
            continue

        if not best:
            print(f"⚠️ Skipping {sample}: No matching primers detected.")
            continue

        fwd_primers = expand_degenerate(best[0])
        rev_primers = expand_degenerate(best[1])
        output_dir = os.path.join(output_base, mode.lower())
        run_cutadapt(fwd, rev, fwd_primers, rev_primers, output_dir)

if __name__ == "__main__":
    base_input = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "uploads"))
    base_output = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "data"))

    for mode in ["16S", "ITS", "AMF", "18S"]:
        if mode in ["AMF", "18S"]:
            folder_name = "18s"
        else:
            folder_name = mode.lower()

        subfolder = os.path.join(base_input, folder_name)

        if not os.path.exists(subfolder):
            print(f"⏩ Skipping {mode}: folder not found at {subfolder}")
            continue

        files = [f for f in os.listdir(subfolder) if f.endswith(".fastq") or f.endswith(".fastq.gz")]
        if not files:
            print(f"⏩ Skipping {mode}: no FASTQ files found in {subfolder}")
            continue

        process_folder(subfolder, mode, output_base=base_output)
