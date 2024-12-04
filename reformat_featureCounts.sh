#!/bin/bash

#SBATCH --job-name="reformatting_featureCounts"       # Job name
#SBATCH --output=featurecounts_%j.out                 # Standard output log
#SBATCH --error=featurecounts_%j.err                  # Standard error log
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1                             # Use 1 CPUs
#SBATCH --time=00:05:00                               # 5 minutes runtime
#SBATCH --mem=1GB                                     # 1GB Memory allocation
#SBATCH --partition=pibu_el8

USER="ecapan"                                         # Define USER variable for general use

# Define paths to input and output directories:
input_file="/data/users/${USER}/rna_seq/RNA-sequencing/featurecounts_output/gene_counts.txt"                   # Where the featureCounts file (.txt) is located.
output_file="/data/users/${USER}/rna_seq/RNA-sequencing/featurecounts_output/reformatted_counts.txt"           # Where the reformatted file (output) will be stored.

tail -n +2 ${input_file} | cut -f 1,7- > ${output_file}

# tail -n +2: Removes the first line.
# cut -f 1,7-: Keeps the first column (gene IDs) and columns 7 onward (sample counts).