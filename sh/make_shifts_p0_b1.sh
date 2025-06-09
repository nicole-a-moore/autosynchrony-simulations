#!/bin/bash
#SBATCH --account=def-jsunday
#SBATCH --nodes=1               
#SBATCH --cpus-per-task=4       
#SBATCH --mem-per-cpu=32GB
#SBATCH --time=40:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user='nicole.moore@mail.mcgill.ca'

# Echo current state
echo "Running job from directory: $(pwd)"
echo "Loading modules..."

module load StdEnv/2020 gcc/9.3.0 udunits/2.2.28  gdal/3.5.1 r/4.2.1 netcdf/4.7.4

echo "Modules loaded."
which R
echo "Checking if R script exists..."
ls -l R/make_shifts.R

# Set R library path
export R_LIBS=~/local/R_libs/

echo "Running R script..."
R CMD BATCH --no-save --no-restore R/make_shifts_p0_b1.R  make_shifts_p0_b1.Rout
echo "Finished running R script."
