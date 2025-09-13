#!/bin/bash
#SBATCH --job-name=tumor_sim
#SBATCH --output=output.log
#SBATCH --error=error.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=72
#SBATCH --time=3:30:00
#SBATCH --mem=8G

module load gcc/12 openmpi/4

echo "Starting ArmA simulation at $(date)"
./project config/PhysiCell_settings_ArmB.xml

echo "Finished both simulations at $(date)"





