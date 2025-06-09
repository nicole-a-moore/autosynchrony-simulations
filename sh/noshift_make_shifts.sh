#!/bin/bash
#SBATCH --account=def-jsunday
#SBATCH --nodes=1               
#SBATCH --cpus-per-task=4       
#SBATCH --mem-per-cpu=32GB
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user='nicole.moore@mail.mcgill.ca'

module load StdEnv/2020 gcc/9.3.0 udunits/2.2.28  gdal/3.5.1 r/4.2.1 netcdf/4.7.4
export R_LIBS=~/local/R_libs/
R CMD BATCH --no-save --no-restore R/noshift_make_shifts.R 
