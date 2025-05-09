#!/usr/bin/env python3

import os
import shutil
import subprocess
import time
from datetime import datetime
import pandas as pd

# === CONFIGURATION ===
PROJECT_DIR = os.path.abspath(".")
EXECUTABLE = os.path.join(PROJECT_DIR, "project")
CONFIG_XML = os.path.join(PROJECT_DIR, "config", "PhysiCell_settings.xml")
OUTPUT_DIR = os.path.join(PROJECT_DIR, "output")
RUNS_DIR = os.path.join(PROJECT_DIR, "runs")
LOG_FILE = os.path.join(RUNS_DIR, "run_log.csv")
# =====================

os.makedirs(RUNS_DIR, exist_ok=True)

print("▶️ Running simulation...")
start_time = time.time()
subprocess.run([EXECUTABLE], check=True)
end_time = time.time()
runtime = round(end_time - start_time, 2)
print(f"✅ Simulation finished in {runtime} seconds.")

timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
run_folder = os.path.join(RUNS_DIR, timestamp)
os.makedirs(run_folder)

# Archive output
for file in os.listdir(OUTPUT_DIR):
    shutil.move(os.path.join(OUTPUT_DIR, file), os.path.join(run_folder, file))

# Save config used
shutil.copy(CONFIG_XML, os.path.join(run_folder, "PhysiCell_settings.xml"))

# Initialize defaults
total_cells = dead_cells = avg_volume = -1
total_tumor = dead_tumor = total_immune = -1
attack_damage_rate = -1.0

# Load results
cells_csv = os.path.join(run_folder, "cells.csv")
if os.path.exists(cells_csv):
    df = pd.read_csv(cells_csv)
    total_cells = len(df)
    
    # Safe accessors
    dead_cells = df["current_phase"].eq(100).sum() if "current_phase" in df.columns else -1
    avg_volume = df["total_volume"].mean() if "total_volume" in df.columns else -1

    if "type" in df.columns:
        total_tumor = df["type"].eq(0).sum()  # type 0 = tumor
        total_immune = df["type"].eq(1).sum()  # type 1 = immune
        dead_tumor = df[(df["type"] == 0) & (df.get("current_phase", pd.Series(-1)).eq(100))].shape[0]

# Read attack_damage_rate from config
try:
    import xml.etree.ElementTree as ET
    tree = ET.parse(CONFIG_XML)
    root = tree.getroot()
    for param in root.findall(".//user_parameters//parameter"):
        if param.get("name") == "attack_damage_rate":
            attack_damage_rate = float(param.get("value"))
except Exception as e:
    print(f"⚠️ Could not parse attack_damage_rate: {e}")

# === Log this run ===
entry = {
    "timestamp": timestamp,
    "total_cells": total_cells,
    "dead_cells": dead_cells,
    "avg_volume": round(avg_volume, 2) if avg_volume != -1 else -1,
    "runtime_s": runtime,
    "attack_damage_rate": attack_damage_rate,
    "total_tumor": total_tumor,
    "dead_tumor": dead_tumor,
    "total_immune": total_immune
}

if os.path.exists(LOG_FILE):
    log_df = pd.read_csv(LOG_FILE)
    log_df = pd.concat([log_df, pd.DataFrame([entry])], ignore_index=True)
else:
    log_df = pd.DataFrame([entry])

log_df.to_csv(LOG_FILE, index=False)
print(f"📄 Logged run to {LOG_FILE}")
