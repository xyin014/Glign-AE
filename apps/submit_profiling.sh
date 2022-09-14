#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=430G
#SBATCH --time=7-00:00:00
#SBATCH --job-name="ae-profile"
#SBATCH -p batch # This is the default partition, you can use any of the following; intel, batch, highmem, gpu

# Print current date
date
hostname
module load centos
centos.sh /rhome/xyin014/Repo/Glign-AE/apps/run_profilings.sh

