#!/bin/bash

#SBATCH --job-name="featurecounts"       # Job name
#SBATCH --output=featurecounts_%j.out    # Standard output log
#SBATCH --error=featurecounts_%j.err     # Standard error log
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4                # Use 4 CPUs
#SBATCH --time=15:00:00                  # 15 hours runtime
#SBATCH --mem=10GB                       # Memory allocation
#SBATCH --partition=pibu_el8

USER="ecapan"                            # Define USER variable for general use

# Define paths to input and output directories:
input_dir="/data/users/${USER}/rna_seq/RNA-sequencing/sorted_mapped_bam"                        # Where the .bam files are located.
output_dir="/data/users/${USER}/rna_seq/RNA-sequencing/featurecounts_output"                    # Where the output will be stored.
annotation_file="/data/users/${USER}/rna_seq/RNA-sequencing/req_mapping_files/homo_sapiens.gtf" # Path of the .gtf file.

# Create the output directory if it doesn't already exist:
mkdir -p ${output_dir}

# List all sorted BAM files in the input directory and store them in an array:
bam_files=(${input_dir}/*_sorted.bam)      # Store BAM file paths in an array.

# Define the output file name:
counts_output="${output_dir}/gene_counts.txt"

# Run featureCounts:
apptainer exec --bind /data/ /containers/apptainer/subread_2.0.1-hed695b0_0.sif \
  featureCounts -T 4 -t exon -g gene_id \
  -a ${annotation_file} \
  -o ${counts_output} \
  ${bam_files[@]}

# -T 4: Use 4 threads (matches SLURM settings).
# -t exon: Count reads mapped to exon features.
# -g gene_id: Group reads by gene ID (as specified in the GTF file's "gene_id" attribute).
# -a: Path to the reference annotation file (GTF format).
# -o: Output file path for gene counts.
# ${bam_files}: All sorted BAM files as input.