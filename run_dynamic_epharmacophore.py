"""
#########################################################################################################################
Title: Automated Dynamic e-Pharmacophore Generation Pipeline

Developed by: Mine Isaoglu, Ph.D.  
Principal Investigator: Serdar Durdagi, Ph.D.  
Affiliation: Computational Drug Design Center (HITMER), Faculty of Pharmacy, Bahçeşehir University, Istanbul, Turkey

Version: May 2025
#########################################################################################################################
"""

import os
import subprocess
import sys
import glob
import shutil
import csv
import time
import argparse
from math import ceil
from multiprocessing import Pool
from schrodinger.structure import StructureReader

# ===========================
# USER ARGUMENTS
# ===========================
def get_default_cores():
    total_cores = os.cpu_count() or 4
    return max(1, int(total_cores * 0.75))  # use 75% of system cores safely

parser = argparse.ArgumentParser(description="Dynamic e-Pharmacophore Workflow from MD frames")
parser.add_argument("--start", type=int, required=True, help="Start index of .mae file (inclusive)")
parser.add_argument("--end", type=int, required=True, help="End index of .mae file (inclusive)")
parser.add_argument("--step", type=int, default=1, help="Step size between frames")
parser.add_argument("--ncores", type=int, default=get_default_cores(), help="Number of CPU cores to use for e-Pharmacophore generation")
parser.add_argument("--batch", type=int, default=get_default_cores(), help="Number of files to process per batch")
args = parser.parse_args()

START_IDX = args.start
END_IDX = args.end
STEP = args.step
USER_NCORES = args.ncores
BATCH_SIZE = args.batch

# ===========================
# PATH INITIALIZATION
# ===========================
BASE_DIR = os.getcwd()
ANALYSIS_DIR = os.path.join(BASE_DIR, "DYNOPHORE_ANALYSIS")
INPUT_DIR = os.path.join(BASE_DIR, "input_mae_files")
PROCESSED_DIR = os.path.join(ANALYSIS_DIR, "PROCESSED_FILES")
HYPOTHESIS_DIR = os.path.join(ANALYSIS_DIR, "saved_HYPOTHESIS")

os.makedirs(ANALYSIS_DIR, exist_ok=True)
os.makedirs(PROCESSED_DIR, exist_ok=True)
os.makedirs(HYPOTHESIS_DIR, exist_ok=True)

SCHRODINGER = os.environ.get("SCHRODINGER18") or os.environ.get("SCHRODINGER18_4") or "/opt/schrodinger2018-4"
if not os.path.isdir(SCHRODINGER):
    print("Schrödinger path not found. Exiting.")
    sys.exit(1)

prepwizard = os.path.join(SCHRODINGER, "utilities", "prepwizard")
pv_convert = os.path.join(SCHRODINGER, "run")
pv_script = "pv_convert.py"
glide_gridgen = os.path.join(SCHRODINGER, "utilities", "generate_glide_grids")
epharmacophores = os.path.join(SCHRODINGER, "utilities", "epharmacophores")

# ===========================
# LOAD & FILTER INPUT FILES
# ===========================
all_maes = sorted(glob.glob(os.path.join(INPUT_DIR, "*.mae")), key=lambda x: int(os.path.splitext(os.path.basename(x))[0]))
selected_maes = [f for f in all_maes if START_IDX <= int(os.path.splitext(os.path.basename(f))[0]) <= END_IDX and (int(os.path.splitext(os.path.basename(f))[0]) - START_IDX) % STEP == 0]

if not selected_maes:
    print("No .mae files found within specified range. Exiting.")
    sys.exit(1)

def nice_run(cmd):
    return subprocess.run(["nice", "-n", "10"] + cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

def process_frame(mae_path):
    base_index = os.path.splitext(os.path.basename(mae_path))[0]
    print(f"[{base_index}] Started processing frame.")

    os.chdir(PROCESSED_DIR)
    work_mae = os.path.join(PROCESSED_DIR, f"{base_index}.mae")
    shutil.copy(mae_path, work_mae)

    print(f"[{base_index}] Running PrepWizard...")
    temp_out = f"{base_index}_pH7.4_prepared.mae"
    final_out = f"{base_index}_prepared.mae"
    try:
        nice_run([prepwizard, f"{base_index}.mae", temp_out,
                  "-NOJOBID", "-noimpref", "-noepik",
                  "-propka_pH", "7.4", "-keepfarwat"])
        os.rename(temp_out, final_out)
        print(f"[{base_index}] PrepWizard completed.")
        subprocess.run(f"rm -f {base_index}-protassign.log {base_index}-protassign.mae {base_index}-protassign-out.mae", shell=True)
    except subprocess.CalledProcessError as e:
        print(f"[{base_index}] PrepWizard failed: {e}")
        return

    print(f"[{base_index}] Splitting ligand and receptor...")
    lig = f"{base_index}_prepared-out_lig.mae"
    rec = f"{base_index}_prepared-out_recep.mae"
    try:
        nice_run([pv_convert, pv_script, "-mode", "split_ligand", final_out, "-o", lig])
        nice_run([pv_convert, pv_script, "-mode", "split_receptor", final_out, "-o", rec])
    except subprocess.CalledProcessError as e:
        print(f"[{base_index}] Splitting failed: {e}")
        return

    print(f"[{base_index}] Calculating ligand center...")
    try:
        with StructureReader(lig) as reader:
            struct = next(reader)
            xs = [a.x for a in struct.atom]
            ys = [a.y for a in struct.atom]
            zs = [a.z for a in struct.atom]
            center = f"{sum(xs)/len(xs):.2f},{sum(ys)/len(ys):.2f},{sum(zs)/len(zs):.2f}"
        print(f"[{base_index}] Ligand center: {center}")
    except Exception as e:
        print(f"[{base_index}] Failed to calculate center: {e}")
        return

    csv_file = f"{base_index}_grid_input.csv"
    with open(csv_file, "w", newline='') as f:
        writer = csv.DictWriter(f, fieldnames=["rec_file", "cent_coor", "hbond_cons", "lig_asl", "res_asl"])
        writer.writeheader()
        writer.writerow({"rec_file": rec, "cent_coor": center, "hbond_cons": "", "lig_asl": "", "res_asl": ""})

    print(f"[{base_index}] Generating Glide grid...")
    try:
        nice_run([glide_gridgen, csv_file])
        print(f"[{base_index}] Glide grid generation complete.")
    except subprocess.CalledProcessError as e:
        print(f"[{base_index}] Glide grid generation failed: {e}")
        return

    actual_zip = None
    while not actual_zip:
        zips = sorted(glob.glob("*.zip"), key=os.path.getmtime, reverse=True)
        for z in zips:
            if "gridgen" in z:
                actual_zip = z
                break
        time.sleep(2)
    gridzip = f"{base_index}_grid.zip"
    shutil.copy(actual_zip, gridzip)

    hypo_prefix = f"{base_index}_hypo"
    eph_input = ["-rec_file", rec, "-lig_file", lig, f"-site_center={center}", "-in_place"]

    print(f"[{base_index}] Generating e-Pharmacophore hypothesis...")
    try:
        nice_run([
            epharmacophores, "-WAIT",
            *eph_input,
            "-fd", "",
            "-f", "7",
            "-site_dist", "2.0",
            "-pair_dist", "4.0",
            "-xvol", "-scale", "0.5",
            "-buff", "2.0", "-limit", "5.0",
            "-HOST", "localhost:1",
            "-j", hypo_prefix
        ])
        print(f"[{base_index}] Hypothesis generated.")
    except subprocess.CalledProcessError as e:
        print(f"[{base_index}] e-Pharmacophore generation failed: {e}")
        return

# ===========================
# BATCH PROCESSING 
# ===========================
total_batches = ceil(len(selected_maes) / BATCH_SIZE)

for batch_idx, batch_start in enumerate(range(0, len(selected_maes), BATCH_SIZE)):
    batch = selected_maes[batch_start:batch_start + BATCH_SIZE]
    print(f"Processing batch {batch_idx + 1}/{total_batches}")

    pool = Pool(processes=min(max(1, USER_NCORES // 2), len(batch)))
    pool.map(process_frame, batch)
    pool.close()
    pool.join()

    print(f"Batch {batch_idx + 1}/{total_batches} completed.")

    # Clean up finished Schrödinger jobs
    try:
        subprocess.run([os.path.join(SCHRODINGER, "jobcontrol"), "-delete", "finished"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        print("[JobControl] Finished jobs cleaned up.")
    except Exception as e:
        print(f"[JobControl] Cleanup failed: {e}")

# ===========================
# FINAL STEP - Collect Hypotheses
# ===========================
for file in glob.glob(os.path.join(PROCESSED_DIR, "*.phypo")):
    shutil.copy(file, os.path.join(HYPOTHESIS_DIR, os.path.basename(file)))

print("All batch analyses completed.")
