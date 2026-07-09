from pathlib import Path
import subprocess

import biom
import h5py


def relabel_fasta(input_fasta: Path, output_fasta: Path, mapping_file: Path, prefix: str) -> dict:
    """Relabel FASTA headers to PREFIX_ASV# and write a mapping file (new_id -> old_id)."""
    asv_map = {}
    counter = 1
    current_seq = ""
    old_id = None

    with input_fasta.open("r", encoding="utf-8", errors="ignore") as infile, \
         output_fasta.open("w", encoding="utf-8") as outfile, \
         mapping_file.open("w", encoding="utf-8") as mapfile:

        for line_num, line in enumerate(infile, 1):
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                # flush previous record
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

        # flush last record
        if old_id and current_seq:
            new_id = f"{prefix}_ASV{counter}"
            outfile.write(f">{new_id}\n{current_seq}\n")
            mapfile.write(f"{new_id}\t{old_id}\n")
            asv_map[old_id] = new_id

    print(f"✅ Cleaned FASTA written to {output_fasta}")
    return asv_map


def relabel_biom(input_biom_fp: Path, output_biom_fp: Path, id_map: dict) -> None:
    """Relabel observation IDs in BIOM using BIOM's update_ids (safer than rebuilding)."""
    print(f"🔁 Relabeling BIOM table: {input_biom_fp.name}")
    table = biom.load_table(str(input_biom_fp))

    # Only update IDs that exist in the table
    # (update_ids ignores non-matching keys, but keeping it clean helps debugging)
    present = set(table.ids(axis="observation"))
    filtered_map = {k: v for k, v in id_map.items() if k in present}

    table.update_ids(filtered_map, axis="observation", inplace=True)

    with h5py.File(str(output_biom_fp), "w") as f:
        table.to_hdf5(f, "Cleaned table with ASV IDs")

    print(f"✅ Saved cleaned BIOM table (HDF5): {output_biom_fp.name}")


def convert_biom_to_tsv(biom_fp: Path, tsv_fp: Path) -> None:
    """Convert biom file to TSV using biom CLI."""
    print(f"📤 Converting {biom_fp.name} to TSV...")
    cmd = [
        "biom", "convert",
        "--input-fp", str(biom_fp),
        "--output-fp", str(tsv_fp),
        "--to-tsv"
    ]
    subprocess.run(cmd, check=True)
    print(f"✅ TSV written to: {tsv_fp}")


if __name__ == "__main__":
    # Add 18s (you export it now)
    datasets = ["16s", "its", "its_amf", "18s", "18s_amf"]

    for domain in datasets:
        path = Path(f"../data/exported/{domain}")
        print(f"\n📦 Cleaning ASVs for: {domain.upper()}")

        if not path.exists():
            print(f"⚠️ Skipping {domain.upper()}: export folder not found: {path}")
            continue

        fasta_in = path / "dna-sequences.fasta"
        biom_in = path / "feature-table.biom"

        fasta_out = path / "cleaned-dna-sequences.fasta"
        mapping_tsv = path / "asv_mapping.tsv"
        biom_out = path / "cleaned-feature-table.biom"
        biom_tsv = path / "cleaned-feature-table.tsv"

        # Require FASTA for mapping (otherwise nothing to relabel)
        if not fasta_in.exists():
            print(f"⚠️ Skipping {domain.upper()}: {fasta_in.name} not found.")
            continue

        # Step 1: Relabel FASTA (creates the ID map)
        asv_map = relabel_fasta(fasta_in, fasta_out, mapping_tsv, prefix=domain.upper())

        # Step 2: Relabel BIOM and convert to TSV (if BIOM exists)
        if biom_in.exists():
            relabel_biom(biom_in, biom_out, id_map=asv_map)
            convert_biom_to_tsv(biom_out, biom_tsv)
        else:
            print(f"⚠️ Skipping BIOM relabel for {domain.upper()}: {biom_in.name} not found.")
