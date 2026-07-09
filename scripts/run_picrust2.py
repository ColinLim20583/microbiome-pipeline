import logging
from pathlib import Path
import pandas as pd
import time
import seaborn as sns
import matplotlib.pyplot as plt
from bioblend.galaxy import GalaxyInstance
import requests

# ------------------- CONFIG -------------------
API_KEY = 'e8abe4f3b5218ac05d555a57ff79d09f'
GALAXY_URL = 'https://usegalaxy.eu'
DATASET_FOLDERS = ["16s", "its", "its_amf", "18s_amf"]
PLOT_OUTPUTS = True  # Toggle plotting for performance

REFERENCE_MAP = {
    "16s": "Prokaryotic 16S rRNA gene",
    "its": "Fungal ITS (only epa-ng)",
    "its_amf": "Fungal 18S (only for epa-ng)",
    "18s_amf": "Fungal 18S (only for epa-ng)"
}

# ------------------- SETUP -------------------
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

# Resolve project directory
if '__file__' in globals():
    PROJECT_DIR = Path(__file__).resolve().parent.parent
else:
    PROJECT_DIR = Path.cwd()

EXPORT_BASE_DIR = PROJECT_DIR / "data" / "exported"
gi = GalaxyInstance(url=GALAXY_URL, key=API_KEY)

# ------------------- TOOL DISCOVERY -------------------
logging.info("Searching for PICRUSt2 Full Pipeline tool...")
tool_id = None
for tool in gi.tools.get_tools():
    if "picrust2 full pipeline" in tool['name'].lower():
        tool_id = tool['id']
        logging.info(f"Found tool: {tool['name']}")
        break

if not tool_id:
    raise RuntimeError("PICRUSt2 Full Pipeline tool not found.")

logging.info(f"Connected as: {gi.users.get_current_user()['email']}")

# ------------------- FILE AND JOB UTILS -------------------
def wait_for_dataset(dataset_id, timeout=900, interval=5):
    start = time.time()
    while time.time() - start < timeout:
        state = gi.datasets.show_dataset(dataset_id)['state']
        if state == 'ok':
            return
        elif state == 'error':
            raise RuntimeError(f"Upload failed for dataset {dataset_id}")
        time.sleep(interval)
    raise TimeoutError(f"Timeout while waiting for dataset {dataset_id} to finish.")

def wait_for_job(job_id, timeout=3600, interval=10):
    start = time.time()
    while time.time() - start < timeout:
        try:
            state = gi.jobs.show_job(job_id)['state']
        except requests.RequestException as e:
            logging.warning(f"Connection issue: {e}, retrying in {interval}s...")
            time.sleep(interval)
            continue
        logging.info(f"Job status: {state}")
        if state in ['ok', 'error']:
            return state
        time.sleep(interval)
    raise TimeoutError(f"Timeout while waiting for job {job_id}")

# ------------------- MAIN LOOP -------------------
for dataset in DATASET_FOLDERS:
    logging.info(f"Processing dataset: {dataset}")
    folder = EXPORT_BASE_DIR / dataset
    fasta_fp = folder / "cleaned-dna-sequences.fasta"
    biom_fp = folder / "cleaned-feature-table.biom"
    map_fp = folder / "asv_mapping.tsv"

    if not (fasta_fp.exists() and biom_fp.exists() and map_fp.exists()):
        logging.warning(f"Skipping {dataset} due to missing files.")
        continue

    reference_data = REFERENCE_MAP.get(dataset)
    if not reference_data:
        logging.warning(f"No reference mapping found for {dataset}. Skipping.")
        continue

    asv_map = pd.read_csv(map_fp, sep='\t', header=None, names=["ASV", "Sequence_ID"])
    asv_dict = dict(zip(asv_map['Sequence_ID'], asv_map['ASV']))

    history = gi.histories.create_history(name=f"PICRUSt2 - {dataset}")
    history_id = history['id']
    logging.info(f"History created: {history['name']}")

    logging.info("Uploading files...")
    try:
        fasta_upload = gi.tools.upload_file(str(fasta_fp), history_id)
        biom_upload = gi.tools.upload_file(str(biom_fp), history_id)
        fasta_id = fasta_upload['outputs'][0]['id']
        biom_id = biom_upload['outputs'][0]['id']
    except Exception as e:
        logging.error(f"File upload failed: {e}")
        continue

    try:
        wait_for_dataset(fasta_id)
        wait_for_dataset(biom_id)
    except Exception as e:
        logging.error(f"Dataset error: {e}")
        continue

    # Input config
    common_inputs = {
        'study_sequences': {'src': 'hda', 'id': fasta_id},
        'abundance_table': {'src': 'hda', 'id': biom_id},
        'placement_tool': 'epa-ng',
        'reference_data': reference_data,
        'min_align': 0.8,
        'trait_tables': 'Default trait table(s)',
        'trait_table_precalc': ['EC', 'KO'],
        'hsp_method': 'mp',
        'transition_cost_weight': 0.5,
        'skip_nsti': 'false',
        'input_type': 'ASV',
        'min_reads': '1',
        'min_samples': '1',
        'stratified': 'true',
        'skip_norm': 'false',
        'reaction_func': 'EC'
    }

    if dataset == "16s":
        common_inputs.update({
            'no_pathways': 'false',
            'skip_minpath': 'false',
            'no_gap_fill': 'false',
            'no_regroup': 'false'
        })
    else:
        common_inputs.update({
            'no_pathways': 'true',
            'skip_minpath': 'true',
            'no_gap_fill': 'true',
            'no_regroup': 'true'
        })

    logging.info("Submitting Galaxy job...")
    try:
        job = gi.tools.run_tool(history_id, tool_id, common_inputs)
        if 'jobs' not in job or not job['jobs']:
            raise RuntimeError("Job submission failed.")
        job_id = job['jobs'][0]['id']
    except Exception as e:
        logging.error(f"Job submission error: {e}")
        continue

    logging.info("Waiting for job to finish...")
    try:
        job_status = wait_for_job(job_id)
        if job_status != 'ok':
            logging.error(f"Job failed for {dataset}")
            continue
    except Exception as e:
        logging.error(f"Job timeout or error: {e}")
        continue

    logging.info("Downloading outputs...")
    try:
        outputs = gi.histories.show_history(history_id, contents=True)
        stratified_fp = None

        for output in outputs:
            out_info = gi.datasets.show_dataset(output['id'])
            if out_info['state'] != 'ok':
                continue
            safe_name = out_info['name'].replace(" ", "_").replace(":", "").replace("/", "-")
            out_path = folder / f"{safe_name}.tsv"
            gi.datasets.download_dataset(out_info['id'], file_path=str(out_path), use_default_filename=False)
            logging.info(f"Saved: {out_path.name}")
            if "stratified_metagenome" in safe_name.lower():
                stratified_fp = out_path
    except Exception as e:
        logging.error(f"Failed to download or parse outputs: {e}")
        continue

    # Analyze stratified output
    if stratified_fp and stratified_fp.exists():
        try:
            df = pd.read_csv(stratified_fp, sep='\t', comment="#")
            df.columns = ["function", "taxonomy"] + list(df.columns[2:])
            df['taxonomy'] = df['taxonomy'].apply(lambda x: asv_dict.get(x, f"[unmapped:{x}]"))
            relabeled_fp = stratified_fp.with_name("Relabeled_" + stratified_fp.name)
            df.to_csv(relabeled_fp, sep='\t', index=False)
            logging.info(f"Relabeled output saved: {relabeled_fp.name}")

            melted = df.melt(id_vars=["function", "taxonomy"], var_name="sample", value_name="abundance")
            melted = melted[melted["taxonomy"] != "[unmapped:unclassified]"]

            summary = (
                melted
                .groupby(["function", "taxonomy"])["abundance"]
                .sum()
                .reset_index()
                .sort_values("abundance", ascending=False)
                .head(20)
            )
            summary_fp = folder / "Top_taxa_function_summary.tsv"
            summary.to_csv(summary_fp, sep='\t', index=False)
            logging.info(f"Summary saved: {summary_fp.name}")

            if PLOT_OUTPUTS:
                summary["label"] = summary["taxonomy"] + " → " + summary["function"]
                plt.figure(figsize=(10, 6))
                sns.barplot(data=summary, y="label", x="abundance", palette="viridis")
                plt.xlabel("Total Predicted Abundance")
                plt.ylabel("Taxon → Function")
                plt.title(f"Top 20 Taxon-Function Contributions ({dataset})")
                plt.tight_layout()
                plot_fp = folder / "Top_taxa_function_summary_plot.png"
                plt.savefig(plot_fp, dpi=300)
                plt.close()
                logging.info(f"Plot saved: {plot_fp.name}")

        except Exception as e:
            logging.warning(f"Error analyzing stratified file: {e}")
    else:
        logging.warning("Stratified output not found.")

logging.info("✅ All datasets processed.")
