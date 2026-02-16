#!/bin/bash
#SBATCH --account=def-jsunday
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --array=1-120
#SBATCH --output=logs/%A_%a.out
#SBATCH --error=logs/%A_%a.err
#SBATCH --time=300:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user='nicole.moore@mail.mcgill.ca'

# Echo current state
echo "Running job from directory: $(pwd)"
echo "Loading modules..."

module load StdEnv/2020 gcc/9.3.0 udunits/2.2.28  gdal/3.5.1 r/4.2.1 netcdf/4.7.4

echo "Modules loaded."
which R
echo "Checking if R script exists..."
ls -l run_one_sim.R

echo "Running R script..."
R CMD BATCH --no-save --no-restore run_one_sim.R 
echo "Finished running R script."
