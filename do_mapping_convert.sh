#!/bin/bash

#SBATCH --job-name="hisat2_mapping_to_bam"       # Job name
#SBATCH --output=hisat2_mapping_%A_%a.out        # Standard output log with job array index
#SBATCH --error=hisat2_mapping_%A_%a.err         # Standard error log with job array index
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8                        # Number of CPUs per task is 8 (same as for HISAT2)
#SBATCH --time=14:00:00                          # 14 hours runtime
#SBATCH --mem=18GB                               # 18GB memory
#SBATCH --partition=pibu_el8
#SBATCH --array=0-11                             # Array range for 12 pairs of fastq files (0-11 for 12 files)

USER="ecapan"                                    # Define variable USER for general use

# Define the input directory where the reads are stored:
input_dir="/data/courses/rnaseq_course/breastcancer_de/reads"

# Define the output directory where the BAM files will be stored:
output_dir="/data/users/${USER}/rna_seq/RNA-sequencing/unsorted_mapped_bam"

# Create the output directory if it doesn't already exist:
mkdir -p ${output_dir}

# Define the genome index:
genome_index="/data/users/${USER}/rna_seq/RNA-sequencing/req_mapping_files/genome_index/genome_index"
# The second "genome_index" (at the end) stands for the prefix of the files obtained from HISAT2 indexing.
# So, the whole path is the path to the index prefix.

# Go to the input directory:
cd ${input_dir}

# Get list of R1 files (files ending in _R1.fastq.gz):
R1_files=($(ls *R1.fastq.gz))

# Get the corresponding R2 files by replacing _R1 with _R2.
# The @ ensures that all elements in the array are processed individually:
for R1_file in "${R1_files[@]}"; do
  R2_file="${R1_file/R1/R2}"

  # Run HISAT2 and pipe the output directly to Samtools to convert SAM to BAM inside the Apptainer container:
  apptainer exec --bind /data/ /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif \
    hisat2 -p 8 -x ${genome_index} -1 ${input_dir}/${R1_file} -2 ${input_dir}/${R2_file} | \
    samtools view -hbS - > ${output_dir}/${R1_file%_R1.fastq.gz}_mappedReads.bam
  # -p 8 tells HISAT2 to use 8 threads (or CPU cores)
  # -x specifies the genome index to use
  # -1 and -2 specify the paired-end read files (R1 and R2)
  # -h displays help (this is optional)
  # -b tells samtools to convert the input into BAM format
  # -S indicates that the input file is in SAM format
  # - (hyphen) signifies that samtools should read from standard input (stdin).
  #   In this case, the samtools command is receiving its input from the pipe (|).
done