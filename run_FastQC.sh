#!/bin/bash

#SBATCH --job-name="JOB_FastQC"        # Job name
#SBATCH --output=fastqc_%j.out         # Standard output log
#SBATCH --error=fastqc_%j.err          # Standard input log
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --time=1:00:00                 # 1h as max. runtime
#SBATCH --mem=2000MB                   # 1000MB memory allocation
#SBATCH --partition=pibu_el8
#SBATCH --array=0-23                   # 24fastq files in total (job 0-23)

# Load modules (in this case apptainer):
module load apptainer

# Define the input and output directories:
input_dir="/data/courses/rnaseq_course/breastcancer_de/reads"
output_dir="/data/users/ecapan/rna_seq/RNA-sequencing/QC_output"

# Create the output directory if it doesn't already exist:
mkdir -p ${output_dir}

# Make an array of input files and select the file for this job:
fastq_files=(${input_dir}/*.fastq.gz)
input_file=${fastq_files[$SLURM_ARRAY_TASK_ID]} # Go through all 24 files.

# Run FastQC on all the FASTQ files:
apptainer exec /containers/apptainer/fastqc-0.12.1.sif fastqc -o ${output_dir} -t 2 ${input_file}
    # CPU threads per task is 2 (-t 2)