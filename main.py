import os
import time
import subprocess
from flask import Flask, request, render_template, redirect, url_for, flash, session
import re

def natural_sort_key(s):
    return [int(text) if text.isdigit() else text.lower() for text in re.split(r'(\d+)', s)]

app = Flask(__name__)
app.secret_key = "supersecret"

BASE_DIR = os.path.abspath(os.path.dirname(__file__))
UPLOAD_FOLDER = os.path.join(BASE_DIR, "uploads")
DATA_FOLDER = os.path.join(BASE_DIR, "data")

MODES = ["16s", "its", "18s"]
SCRIPT_ORDER = [
    "scripts/fastqc.py",            # Step 1: Pre-trim QC
    "scripts/priming.py",           # Step 2: Primer trimming
    "scripts/fastqc.py",            # Step 3: Post-trim QC
    "scripts/generate_metadata.py", # Step 4: Combine metadata
    "scripts/generate_manifest.py",  # ✅ Step 5: Generate manifest for QIIME2
    "scripts/import_qiime2.py",  # ✅ Step 6: Import demux.qza and qzv for pair end
    "scripts/dada2_denoise.py",  # ✅ Step 7: Denoising with dada2 in Qiime2
    "scripts/taxonomic_classification.py",  # Step 8: Run classifiers (16S, ITS, ITS-AMF, 18S-AMF)
    "scripts/export_all_qza.py",  # , ✅ Step 9: create biom table for 16s, ITS, AMF and 18s
    "scripts/clean_asv_fasta.py",  # ✅ Step 10: standardising ASV ID with mapping of feature ID
    "scripts/filter_amf_table.py", # ✅ Step 11: Filter table for AMF (Glomeromycetes)
    "scripts/run_qiime_diversity.py",  # ✅ Step 12: Generating Alpha and beta diversity with shanon index
    "scripts/generate_asv_table.py",  # ✅ Step 13: Generating ASV Table
    "scripts/heatmap.py",  # ✅ Step 14: Generating Heatmap
    "scripts/krona_chart.py",  # ✅ Step 15: Generating Krona chart
    "scripts/export_qiime2.py",  # ✅ Step 16:Exporting qiime2 data into picrust2 format. Then change env to picrust2
    "scripts/generate_functional_order.py"  # ✅ Step 17: Generate order for TND
]

@app.route('/')
def home():
    session.clear()
    return redirect(url_for("upload", mode="16s"))

@app.route("/upload/<mode>", methods=["GET", "POST"])
def upload(mode):
    if mode not in MODES:
        return f"❌ Invalid mode: {mode}", 400

    folder_path = os.path.join(UPLOAD_FOLDER, mode)
    os.makedirs(folder_path, exist_ok=True)

    if request.method == "POST":
        # Handle uploaded FASTQ files
        files = request.files.getlist("fastq_files")
        if files:
            sample_ids = set()
            for f in files:
                if f and f.filename:
                    filename = f.filename
                    f.save(os.path.join(folder_path, filename))
                    name = re.sub(r'\.(fastq|fq)(\.gz)?$', '', filename, flags=re.IGNORECASE)
                    name = re.sub(r'(_R?[12])$', '', name, flags=re.IGNORECASE)  # removes _R1/_R2/_1/_2
                    name = re.sub(r'_trimmed$', '', name, flags=re.IGNORECASE)  # removes trailing _trimmed
                    sample_id = name.strip("_")
                    if sample_id:
                        sample_ids.add(sample_id)
            session[f"{mode}_samples"] = sorted(sample_ids, key=natural_sort_key)
            return redirect(url_for("metadata", mode=mode))

    return render_template("index.html", mode=mode)

@app.route('/metadata/<mode>', methods=["GET", "POST"])
def metadata(mode):
    if mode not in MODES:
        return "Invalid mode", 400

    samples = session.get(f"{mode}_samples", [])
    metadata_file = os.path.join(DATA_FOLDER, f"{mode}_metadata.tsv")

    if request.method == "POST":
        total = int(request.form.get("total", 0))
        with open(metadata_file, "w") as f:
            f.write("sample-id\ttreatment\ttimepoint\tlocation\tsoil\tcrop\tph\tnote\n")
            for i in range(1, total + 1):
                sid = request.form.get(f"sample_id_{i}", "")
                treatment = request.form.get(f"treatment_{i}", "")
                timepoint = request.form.get(f"timepoint_{i}", "")
                location = request.form.get(f"location_{i}", "")
                soil = request.form.get(f"soil_{i}", "")
                crop = request.form.get(f"crop_{i}", "")
                ph = request.form.get(f"ph_{i}", "")
                note = request.form.get(f"note_{i}", "")
                f.write(f"{sid}\t{treatment}\t{timepoint}\t{location}\t{soil}\t{crop}\t{ph}\t{note}\n")

        next_index = MODES.index(mode) + 1
        if next_index < len(MODES):
            return redirect(url_for("confirm_next", next_mode=MODES[next_index]))
        else:
            return redirect(url_for("run_pipeline"))

    return render_template("metadata.html", samples=samples, mode=mode)

@app.route("/confirm/<next_mode>", methods=["GET", "POST"])
def confirm_next(next_mode):
    if next_mode not in MODES:
        return "Invalid next mode", 400

    if request.method == "POST":
        if "yes" in request.form:
            return redirect(url_for("upload", mode=next_mode))
        elif "no" in request.form:
            return advance_or_finish(next_mode)

    return render_template("confirm_next.html", next_mode=next_mode)

@app.route("/run_pipeline")
def run_pipeline():
    return redirect(url_for("processing_status", step=0))

@app.route("/processing_status/<int:step>")
def processing_status(step):
    steps = [
        "Running FastQC on raw reads...",
        "Trimming and primer removal...",
        "Running FastQC on trimmed reads...",
        "Generating combined metadata...",
        "Generating manifest files for QIIME2 import...",
        "Importing Demux pair report...",
        "Denosing with Dada2, it will take sometime, please wait...",
        "Running classifiers...",
        "Creating Biom table and DNA FASTA...",
        "Standardising ASV ID... ",
        "Filtering table for AMF ",
        "Generating Alpha and Beta diversity...",
        "Generate ASV Table...",
        "Generate Heatmap...",
        "Generate Krona Chart...",
        "exporting_qiime2 to picrust2 format...",
        "Generating TND..."
    ]

    if step < len(steps):
        try:
            # Show current step
            current_step_text = steps[step]
            # Actually run the script
            run_script(SCRIPT_ORDER[step])
        except Exception as e:
            return f"<h3>❌ Error in step {step + 1}: {e}</h3>"

        # Now proceed to next step immediately
        return redirect(url_for("processing_status", step=step + 1))

    # If all steps are done
    return render_template("status.html", current="✅ All steps complete.", next_step=0, done=True)

def run_script(script_path):
    full_path = os.path.join(BASE_DIR, script_path)
    if not os.path.exists(full_path):
        raise FileNotFoundError(f"{script_path} not found")
    subprocess.run(["python", full_path], check=True)

def advance_or_finish(next_mode):
    if next_mode == MODES[-1]:
        return redirect(url_for("run_pipeline"))
    else:
        next_index = MODES.index(next_mode) + 1
        return redirect(url_for("confirm_next", next_mode=MODES[next_index]))

if __name__ == "__main__":
    app.run(debug=True)
