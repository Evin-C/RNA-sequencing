#!/bin/bash

#SBATCH --job-name="download_genome_and_annotation" # Job name
#SBATCH --output=fastqc_%j.out         # Standard output log
#SBATCH --error=fastqc_%j.err          # Standard input log
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0:30:00                 # 30min as max. runtime
#SBATCH --mem=1000MB                   # 1000MB memory allocation
#SBATCH --partition=pibu_el8

USER="ecapan"                          # define variable USER for a more general use

# Define the output directory (where the file(s) should be downloaded to):
output_dir="/data/users/${USER}/rna_seq/RNA-sequencing/req_mapping_files"

# Go to the output directory, where the reference genome sequence and the associated annotation file should be downloaded to:
cd ${output_dir}

# Download the latest reference genome sequence (.fa) (currently it is release 113) of our species (human):
wget https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# And download the associated annotation (.gtf) (so also release 113):
wget https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz

# After the download, go to the directory of the downloaded files and check that the downloaded files are intact by computing checksums:
# sum Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz and 
# sum Homo_sapiens.GRCh38.113.gtf.gz
# and comparing them to the values in the CHECKSUMS file on the ftp server.