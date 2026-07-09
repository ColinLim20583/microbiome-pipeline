import os
import pandas as pd
import re

BASE_DIR = os.path.abspath(os.path.dirname(__file__))
DATA_DIR = os.path.join(BASE_DIR, "../data")


def natural_sort_key(s):
    return [int(text) if text.isdigit() else text.lower() for text in re.split(r'(\d+)', s)]


def get_fastq_pairs(folder):
    pairs = {}
    for file in os.listdir(folder):
        if not (file.endswith(".fastq") or file.endswith(".fq")):
            continue

        full_path = os.path.abspath(os.path.join(folder, file))

        match = re.match(r"(.+?)_(R1|R2).*\.f(ast)?q$", file, re.IGNORECASE)
        if not match:
            continue

        sample_id, direction = match.groups()[0], match.groups()[1].upper()

        if sample_id not in pairs:
            pairs[sample_id] = {"R1": None, "R2": None}

        if direction == "R1":
            pairs[sample_id]["R1"] = full_path
        elif direction == "R2":
            pairs[sample_id]["R2"] = full_path

    return pairs


def main():
    for mode in os.listdir(DATA_DIR):
        mode_dir = os.path.join(DATA_DIR, mode)
        if not os.path.isdir(mode_dir):
            continue

        print(f"📁 Processing {mode.upper()} files...")

        pairs = get_fastq_pairs(mode_dir)
        rows = []

        for sample in sorted(pairs.keys(), key=natural_sort_key):
            paths = pairs[sample]
            if paths["R1"] and paths["R2"]:
                rows.append({
                    "sample-id": sample,
                    "forward-absolute-filepath": paths["R1"],
                    "reverse-absolute-filepath": paths["R2"]
                })
            else:
                print(f"⚠️ Incomplete pair for sample: {sample} in {mode}")

        if rows:
            df = pd.DataFrame(rows)
            manifest_path = os.path.join(DATA_DIR, f"{mode}_manifest.tsv")
            df.to_csv(manifest_path, sep="\t", index=False)
            print(f"✅ Manifest written: {manifest_path}")
        else:
            print(f"⚠️ No complete sample pairs found for {mode.upper()}.")


if __name__ == "__main__":
    main()
