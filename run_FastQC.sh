#!/bin/bash

#SBATCH --job-name="JOB_FastQC"        # Job name
#SBATCH --output=fastqc_%j.out         # Standard output log
#SBATCH --error=fastqc_%j.err          # Standard input log
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00                 # 1h as max. runtime
#SBATCH --mem=1000MB                   # 1000MB memory allocation
#SBATCH --partition=pibu_el8
# #SBATCH --array=0-23                 # 24fastq files in total (job 0-23)

# Load modules (we need this file: fastqc-0.12.1.sif)
#source /containers/apptainer/fastqc-0.12.1.sif
module load apptainer

# Define the input and output directories
input_dir="/data/courses/rnaseq_course/breastcancer_de/reads"
output_dir="/data/users/ecapan/rna_seq/RNA-sequencing/QC_output"

# Run FastQC on all the FASTQ files
apptainer exec /containers/apptainer/fastqc-0.12.1.sif fastqc -o ${output_dir} -t 1 ${input_dir}/*.fastq.gz
