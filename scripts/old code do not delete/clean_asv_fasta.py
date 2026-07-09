from pathlib import Path
import biom
import h5py
import subprocess

def relabel_fasta(input_fasta, output_fasta, mapping_file, prefix):
    asv_map = {}
    counter = 1
    current_seq = ""
    old_id = None

    with input_fasta.open("r", encoding="utf-8", errors="ignore") as infile, \
         output_fasta.open("w") as outfile, \
         mapping_file.open("w") as mapfile:

        for line_num, line in enumerate(infile, 1):
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if old_id and current_seq:
                    new_id = f"{prefix}_ASV{counter}"
                    outfile.write(f">{new_id}\n{current_seq}\n")
                    mapfile.write(f"{new_id}\t{old_id}\n")
                    asv_map[old_id] = new_id
                    counter += 1
                old_id = line[1:].strip()
                if not old_id:
                    print(f"❌ Empty header at line {line_num}. Skipping.")
                    old_id = None
                current_seq = ""
            else:
                if old_id:
                    current_seq += line
                else:
                    print(f"❌ Sequence without valid header at line {line_num}: {line}")

        # Final sequence
        if old_id and current_seq:
            new_id = f"{prefix}_ASV{counter}"
            outfile.write(f">{new_id}\n{current_seq}\n")
            mapfile.write(f"{new_id}\t{old_id}\n")
            asv_map[old_id] = new_id

    print(f"✅ Cleaned FASTA written to {output_fasta}")
    return asv_map

def relabel_biom(input_biom_fp, output_biom_fp, id_map):
    print(f"🔁 Relabeling BIOM table: {input_biom_fp.name}")
    table = biom.load_table(str(input_biom_fp))

    new_obs_ids = [id_map.get(i, i) for i in table.ids(axis='observation')]
    new_table = biom.Table(
        table.matrix_data,
        new_obs_ids,
        table.ids(axis='sample'),
        observation_metadata=table.metadata(axis='observation')
    )

    # Save as HDF5 (BIOM v2)
    with h5py.File(str(output_biom_fp), "w") as f:
        new_table.to_hdf5(f, "Cleaned table with ASV IDs")

    print(f"✅ Saved cleaned BIOM table (HDF5): {output_biom_fp.name}")

def convert_biom_to_tsv(biom_fp: Path, tsv_fp: Path):
    """Convert biom file to TSV using biom CLI"""
    print(f"📤 Converting {biom_fp.name} to TSV...")
    cmd = [
        "biom", "convert",
        "--input-fp", str(biom_fp),
        "--output-fp", str(tsv_fp),
        "--to-tsv"
    ]
    subprocess.run(cmd, check=True)
    print(f"✅ TSV written to: {tsv_fp}")

# --- Main loop ---
datasets = ["16s", "its", "its_amf", "18s_amf"]
for domain in datasets:
    path = Path(f"../data/exported/{domain}")
    print(f"\n📦 Cleaning ASVs for: {domain.upper()}")

    # File paths
    fasta_in = path / "dna-sequences.fasta"
    fasta_out = path / "cleaned-dna-sequences.fasta"
    mapping_tsv = path / "asv_mapping.tsv"
    biom_in = path / "feature-table.biom"
    biom_out = path / "cleaned-feature-table.biom"
    biom_tsv = path / "cleaned-feature-table.tsv"

    # Step 1: Relabel FASTA
    asv_map = relabel_fasta(fasta_in, fasta_out, mapping_tsv, prefix=domain.upper())

    # Step 2: Relabel BIOM and convert to TSV
    if biom_in.exists():
        relabel_biom(biom_in, biom_out, id_map=asv_map)
        convert_biom_to_tsv(biom_out, biom_tsv)
    else:
        print(f"⚠️ Skipping BIOM relabel for {domain.upper()}: {biom_in.name} not found.")
