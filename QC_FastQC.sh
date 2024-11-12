#!/bin/bash

#SBATCH --job-name="QC_FastQC"
#SBATCH --mail-type=fail               # Notify when job fails
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00                 # Maximum runtime of 2 hours
#SBATCH --mem=1000MB
#SBATCH --partition=pibu_el8
#SBATCH --array=0-23                   # 24fastq files in total (job 0-23)

