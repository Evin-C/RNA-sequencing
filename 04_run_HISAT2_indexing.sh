#!/bin/bash

#SBATCH --job-name="hisat2_indexing"         # Job name
#SBATCH --output=hisat2_indexing_%j.out      # Standard output log
#SBATCH --error=hisat2_indexing_%j.err       # Standard error log
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4                    # Number of CPUs is 4 (same as for HISAT2)
#SBATCH --time=1:00:00                       # 1 hour runtime
#SBATCH --mem=8GB                            # 8GB memory allocation
#SBATCH --partition=pibu_el8

USER="ecapan"                                # Define variable USER for general use

# Define the directory where you want the index files to be stored:
output_dir="/data/users/${USER}/rna_seq/RNA-sequencing/req_mapping_files/genome_index"

# Create the output directory if it doesn't already exist:
mkdir -p ${output_dir}

# Define the directory where the genome file (genome.fa) is stored:
input_dir="/data/users/${USER}/rna_seq/RNA-sequencing/req_mapping_files"

# Go to the input directory:
cd ${input_dir}

# Run HISAT2 to index the genome using Apptainer:
apptainer exec --bind /data/ /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif hisat2-build -p 4 ${input_dir}/genome.fa ${output_dir}/genome_index
# -p 4 tells HISAT2 to use 4 threads (or CPU cores)
# "genome_index" will be the prefix for the output index files (which will be stored in the same named directory)