#!/bin/bash
#SBATCH --job-name=tumor_sim
#SBATCH --output=output.log
#SBATCH --error=error.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00

module load gcc/12 openmpi/4
module load imagemagick ffmpeg

# Navigate to your working directory
cd /kaznay/PhysiCell

# Compile the project (if not already compiled)
make

# 1. Run the simulation
./project config/PhysiCell_settings_ArmB.xml








