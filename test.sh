#!/bin/bash

#SBATCH -D /net/nfs/home/meignin/Documents/RainSat/

# Job name
#SBATCH -J EDA
#SBATCH --error=/net/nfs/home/meignin/Documents/RainSat/dump/slurm-%j.err
#SBATCH --output=/net/nfs/home/meignin/Documents/RainSat/dump/slurm-%j.out
#SBATCH --time=7-00:00:00

# Load the user profile
source /etc/profile

# Clear any previous module environment
module purge


# Activate the conda environment
anaconda_bin_path=/net/nfs/home/meignin/miniconda3/bin
source "${anaconda_bin_path}/activate" rainsat

echo "Script launch"

# Run the Python script
python3 extraction/rain_rate_distrib.py

# Print confirmation
echo "Script executed successfully"

# Deactivate conda environment
conda deactivate

# Exit the script
exit 0
