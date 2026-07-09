import os
import re
import subprocess
import sys

# All primers / paths / options come from the central config so nothing is
# hard-coded here. Falls back to built-in defaults if config.py is unavailable.
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
try:
    import config as cfg
except Exception:
    cfg = None

if cfg is not None:
    DISCARD_UNTRIMMED = cfg.get("cutadapt", {}).get("discard_untrimmed", True)
    ERROR_RATE = str(cfg.get("cutadapt", {}).get("error_rate", 0.12))
    CUTADAPT_TIMES = str(cfg.get("cutadapt", {}).get("times", 1))
    _P = cfg.PRIMERS
    PRIMERS_16S = tuple(_P["16s"])
    PRIMERS_ITS = tuple(_P["its"])
    PRIMERS_18S = tuple(_P["18s"])
    PRIMERS_AMF = tuple(_P["amf"])
    FASTQ_EXTS = cfg.FASTQ_EXTS
else:
    DISCARD_UNTRIMMED = True
    ERROR_RATE = "0.12"
    CUTADAPT_TIMES = "1"
    PRIMERS_16S = ("ACTCCTACGGGAGGCAGCAG", "GGACTACHVGGGTWTCTAAT")          # 338F / 806R
    PRIMERS_ITS = ("CTTGGTCATTTAGAGGAAGTAA", "GCTGCGTTCTTCATCGATGC")        # ITS1F / ITS2
    PRIMERS_18S = ("CGWTAACGAACGAGACCT", "AICCATTCAATCGGTAIT")              # FF390 / FR1 (has I)
    PRIMERS_AMF = ("ATCAACTTTCGATGGTAGGATAGA", "GAACCCAAACACTTTGGTTTCC")    # AML1 / AML2
    FASTQ_EXTS = (".fastq", ".fastq.gz", ".fq", ".fq.gz")

def normalize_primer(p: str) -> str:
    # Cutadapt does not treat inosine "I" as wildcard; use N (any base)
    return p.upper().replace("I", "N")

def find_pairs(folder):
    pairs = {}
    for fn in sorted(os.listdir(folder)):
        if not fn.lower().endswith(tuple(FASTQ_EXTS)):
            continue

        full = os.path.join(folder, fn)

        m = re.match(r"(.+?)(_R?)([12])(\D.*)?\.(f(ast)?q)(\.gz)?$", fn, flags=re.IGNORECASE)
        if m:
            sample = m.group(1)
            which = m.group(3)
            r = "R1" if which == "1" else "R2"
            pairs.setdefault(sample, {"R1": None, "R2": None})[r] = full
            continue

        m2 = re.match(r"(.+?)_(R1|R2)\.(fastq|fq)(\.gz)?$", fn, flags=re.IGNORECASE)
        if m2:
            sample = m2.group(1)
            r = m2.group(2).upper()
            pairs.setdefault(sample, {"R1": None, "R2": None})[r] = full
            continue

        m3 = re.match(r"(.+?)_([12])\.(fastq|fq)(\.gz)?$", fn, flags=re.IGNORECASE)
        if m3:
            sample = m3.group(1)
            r = "R1" if m3.group(2) == "1" else "R2"
            pairs.setdefault(sample, {"R1": None, "R2": None})[r] = full
            continue

    complete = {s: d for s, d in pairs.items() if d["R1"] and d["R2"]}
    for s, d in pairs.items():
        if not (d["R1"] and d["R2"]):
            print(f"⚠️ Incomplete pair for sample: {s} in {folder}")
    return complete

def run_cutadapt(sample_id, r1_file, r2_file, fwd, rev, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    out_f = os.path.join(output_dir, f"{sample_id}_R1_trimmed.fastq")
    out_r = os.path.join(output_dir, f"{sample_id}_R2_trimmed.fastq")

    cmd = [
        "cutadapt",
        "-e", ERROR_RATE,
        "--times", CUTADAPT_TIMES,
        "-g", fwd,
        "-G", rev,
        "-o", out_f,
        "-p", out_r,
        r1_file, r2_file
    ]

    if DISCARD_UNTRIMMED:
        cmd += ["--discard-untrimmed"]

    print(f"✂️ Running Cutadapt: {sample_id}")
    subprocess.run(cmd, check=True)
    print(f"✅ Cutadapt completed for {sample_id}")

def process_folder(folder, mode, output_base, primer_pair):
    print(f"\n📂 Processing folder: {folder} [{mode}]")
    pairs = find_pairs(folder)
    if not pairs:
        print(f"⚠️ No complete pairs found in {folder}")
        return

    fwd = normalize_primer(primer_pair[0])
    rev = normalize_primer(primer_pair[1])

    output_dir = os.path.join(output_base, mode.lower())
    for sample_id, d in pairs.items():
        print(f"🔎 {sample_id}: forced primers, output={output_dir}")
        run_cutadapt(sample_id, d["R1"], d["R2"], fwd, rev, output_dir)

if __name__ == "__main__":
    base_input = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "uploads"))
    base_output = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "data"))

    # IMPORTANT: AMF and 18S must NOT point to the same folder.
    PLAN = [
        ("16S", "16s", PRIMERS_16S),
        ("ITS", "its", PRIMERS_ITS),
        ("18S", "18s", PRIMERS_18S),
        ("AMF", "amf", PRIMERS_AMF),
    ]

    for mode, folder_name, primers in PLAN:
        subfolder = os.path.join(base_input, folder_name)

        if not os.path.exists(subfolder):
            print(f"⏩ Skipping {mode}: folder not found at {subfolder}")
            continue

        files = [f for f in os.listdir(subfolder) if f.lower().endswith(tuple(FASTQ_EXTS))]
        if not files:
            print(f"⏩ Skipping {mode}: no FASTQ files found in {subfolder}")
            continue

        process_folder(subfolder, mode, output_base=base_output, primer_pair=primers)
