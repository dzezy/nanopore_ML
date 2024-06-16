#!/bin/bash
#SBATCH --job-name=nextflow
#SBATCH --output=nextflow.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
#SBATCH --mem-per-cpu=8G

# activate conda environment for nextflow
source /opt/intel/intelpython3/etc/profile.d/conda.sh
conda activate /home/dzhou/miniconda3/envs/nextflow

# cd into directory with .nf file
cd /home/dzhou/nanopore

# run nextflow pipeline script
time nextflow run main.nf -resume -with-report report.html
