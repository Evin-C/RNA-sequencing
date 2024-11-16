#!/bin/bash

#SBATCH --job-name="unzip_genome_annotation" # Job name
#SBATCH --output=unzip_%j.out               # Standard output log
#SBATCH --error=unzip_%j.err                # Standard error log
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1                   # Number of CPUs is 1
#SBATCH --time=0:15:00                      # 15min as max. runtime
#SBATCH --mem=500MB                         # 500MB memory allocation
#SBATCH --partition=pibu_el8

USER="ecapan"                               # Define variable USER for general use

# Define the directory where the downloaded files (still zipped) are stored:
input_dir="/data/users/${USER}/rna_seq/RNA-sequencing/req_mapping_files"

# Go to the directory:
cd ${input_dir}

# Unzip the genome file (.fa.gz):
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Unzip the annotation file (.gtf.gz):
gunzip Homo_sapiens.GRCh38.113.gtf.gz

# Rename the genome file and GTF file to a simpler to use name:
mv Homo_sapiens.GRCh38.dna.primary_assembly.fa genome.fa
mv Homo_sapiens.GRCh38.113.gtf homo_sapiens.gtf