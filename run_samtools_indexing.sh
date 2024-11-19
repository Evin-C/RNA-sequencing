#!/bin/bash

#SBATCH --job-name="samtools_indexing_bam"    # Job name
#SBATCH --output=index_bam_%A_%a.out    # Standard output log with job array index
#SBATCH --error=index_bam_%A_%a.err     # Standard error log with job array index
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1               # Only 1 CPU needed for indexing
#SBATCH --time=1:00:00                  # 1 hour runtime
#SBATCH --mem=6GB                       # 6GB memory allocation
#SBATCH --partition=pibu_el8
#SBATCH --array=0-11                    # Array range for 12 BAM files (job 0-11)

USER="ecapan"                           # Define USER variable for general use

# Define the input directory:
input_dir="/data/users/${USER}/rna_seq/RNA-sequencing/sorted_mapped_bam"

# Make a list of the sorted BAM files:
bam_files=($(ls ${input_dir}/*_sorted.bam))

# Select current BAM file based on the job array index:
bam_file=${bam_files[$SLURM_ARRAY_TASK_ID]}

# Index the sorted BAM file:
apptainer exec --bind /data/ /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif \
  samtools index ${bam_file}

#The resulting .bai files will be in the same folder as the input .bam files.