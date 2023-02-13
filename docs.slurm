#!/bin/bash
#SBATCH --job-name=lrv_cmp
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=cpuq
#SBATCH --mem=100GB
#SBATCH --output=lrv_cmp_%j.log
#SBATCH --error=lrv_cmp_%j.stderr

module load singularity

make docs_singularity
