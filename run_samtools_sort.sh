#!/bin/bash

#SBATCH --job-name="sort_bam_files"    # Job name
#SBATCH --output=sort_bam_%A_%a.out    # Standard output log with job array index
#SBATCH --error=sort_bam_%A_%a.err     # Standard error log with job array index
#SBATCH --nodes=1
#SBATCH --cpus-per-task=6              # Number of CPUs is 6 (same as for samtools)
#SBATCH --time=2:30:00                 # 2.5h as max. runtime
#SBATCH --mem=30GB                     # 30GB memory allocation
#SBATCH --partition=pibu_el8
#SBATCH --array=0-11                   # Array range for 12 files

USER="ecapan"                          # Define USER variable for general use

# Define the input and output directories:
input_dir="/data/users/${USER}/rna_seq/RNA-sequencing/unsorted_mapped_bam"
output_dir="/data/users/${USER}/rna_seq/RNA-sequencing/sorted_mapped_bam"

# Create the output directory if it doesn't already exist:
mkdir -p ${output_dir}

# Make a list of the unsorted BAM files:
bam_files=($(ls ${input_dir}/*.bam))

# Select the current BAM file based on the job array index:
bam_file=${bam_files[$SLURM_ARRAY_TASK_ID]}

# Define the base name for the output file (basename strips directory and suffix from filenames):
base_name=$(basename ${bam_file} .bam)

# Sort the BAM file by genomic coordinates:
apptainer exec --bind /data/ /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif \
  samtools sort -@ 6 -o ${output_dir}/${base_name}_sorted.bam ${bam_file}
# -@ 6 uses 6 threads for sorting.
# -o specifies the output file.

# Check exit status of samtools and log success/failure:
if [ $? -eq 0 ]; then
  echo "$(date): Sorting completed successfully for ${bam_file}" >> ${output_dir}/sort_success.log
else
  echo "$(date): Sorting failed for ${bam_file}" >> ${output_dir}/sort_error.log
fi
# $? holds the exit status of the most recently executed command (in this case, the exit status of samtools sort).
#    0 means the command completed successfully. Non-zero values (e.g., 1, 2) indicate that something went wrong.
#    After samtools sort runs, successes are written to sort_success.log and failures are written to sort_error.log.
#    With this, silent failures can be avoided. If a file wasn’t sorted properly, it will show in the log file.
# ${date} creates timestamps for better tracking.