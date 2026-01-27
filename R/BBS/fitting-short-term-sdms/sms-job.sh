#!/bin/bash
#SBATCH --account=def-jsunday
#SBATCH --mem-per-cpu=32GB
#SBATCH --time=168:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user='nicole.moore@mail.mcgill.ca'

module load StdEnv/2020 gcc/9.3.0 r/4.2.1 gdal/3.5.1

Rscript run-sdms.R