#!/bin/bash
#SBATCH --job-name=tumor_sim
#SBATCH --output=output.log
#SBATCH --error=error.log
#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=12:00:00


module load gcc/12 openmpi/4
module load imagemagick ffmpeg

# Navigate to your working directory
cd /ptmp/akaynak/PhysiCell || { echo "Directory not found"; exit 1; }

# Compile the project (if not already compiled)
make -j8  # use all available cores for faster build

# 1. Run the simulation
./project config/PhysiCell_settings_ArmB.xml

# 2. Convert SVGs to PNGs (only if SVGs exist)
mkdir -p png_frames
for file in output/snapshot*.svg; do
    if [[ -f "$file" ]]; then
        magick "$file" "png_frames/$(basename "$file" .svg).png"
    fi
done

# 3. Generate movie (only if PNGs exist)
if ls png_frames/snapshot*.png 1> /dev/null 2>&1; then
    ffmpeg -framerate 5 -i png_frames/snapshot%08d.png -pix_fmt yuv420p tumor_simulation.mp4
else
    echo "No PNG files found — skipping movie generation."
fi
