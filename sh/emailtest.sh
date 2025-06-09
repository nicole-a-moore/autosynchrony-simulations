#!/bin/bash
#SBATCH --job-name=test-email
#SBATCH --time=00:01:00
#SBATCH --mem=1G
#SBATCH --mail-user=nicole.moore@mail.mcgill.ca
#SBATCH --mail-type=BEGIN,END,FAIL

echo "Testing email notifications..."
sleep 30