#!/bin/bash

#SBATCH --job-name="trimmed_FastQC"    # Job name
#SBATCH --output=fastqc_%j.out         # Standard output log
#SBATCH --error=fastqc_%j.err          # Standard input log
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4              # Use 4 CPUs - as specified for Trimmomatic threads
#SBATCH --time=1:00:00                 # 1h as max. runtime
#SBATCH --mem=2000MB                   # 2000MB memory allocation
#SBATCH --partition=pibu_el8
#SBATCH --array=0-11                   # 24fastq files in total = 12 pairs (*_R1 and *_R2) (job 0-11)

USER="ecapan"                          # define variable USER for a more general use

# Load modules:
module load Java/11.0.18               # needed to run Trimmomatic
module load Trimmomatic/0.39-Java-11   # trims the 

# Define the input and output directories:
input_dir="/data/courses/rnaseq_course/breastcancer_de/reads"
trimmed_dir="/data/users/${USER}/rna_seq/RNA-sequencing/trimmed_reads"
output_dir="/data/users/${USER}/rna_seq/RNA-sequencing/trimQC_output"

# Create the output directory if they don't already exist:
mkdir -p ${trimmed_dir}
mkdir -p ${output_dir}

# Create an array for R1 files and use the index for the corresponding R2 files:
read1_files=(${input_dir}/*_R1.fastq.gz)    # Array of files with _R1 suffix
READ1=${read1_files[$SLURM_ARRAY_TASK_ID]}  # Get the R1 file for this job
READ2=${READ1/_R1.fastq.gz/_R2.fastq.gz}    # Get the R2 filename by replacing "_R1.fastq.gz" with "_R2.fastq.gz"

# Define output filenames based on the sample name:
sample_name=$(basename ${READ1} _R1.fastq.gz)    # Extract the sample name using "basename" as path
                                                 # and "_R1.fastq.gz" as the suffix to be removed from the end of the file name.
trimmed_file1="${trimmed_dir}/${sample_name}_R1trim.fastq.gz"
trimmed_file2="${trimmed_dir}/${sample_name}_R2trim.fastq.gz"
# Need unpaired files in case of trimming leading to loss of a read in a pair due to poor quality bases:
unpaired_file1="${trimmed_dir}/${sample_name}_R1unpaired.fastq.gz"
unpaired_file2="${trimmed_dir}/${sample_name}_R2unpaired.fastq.gz"

# Now the reads can be trimmed by running Trimmomatic:
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -phred33 -threads 4 \
    ${READ1} ${READ2} \
    ${trimmed_file1} ${unpaired_file1} \
    ${trimmed_file2} ${unpaired_file2} \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50
# "java -jar" runs a Java program from a JAR file (Trimmomatic is implemented as a Java program).
# "$EBROOTTRIMMOMATIC/trimmomatic-0.39.jar" points to the location of the Trimmomatic JAR file.
# "PE" stands for Paired-End trimming, meaning we provide two FASTQ files (R1 and R2) and they
#      will be handled together by Trimmomatic.
# "-phred33" indicates the phred score encoding that will be used in the quality score of the sequences.
#            In "phred33" quality scores range from 0 to 41, with 0 being the worst, and 41 being the best quality.
# "${READ1} ${READ2}" are the input files.
# "${trimmed_file1} ${unpaired_file1} ${trimmed_file2} ${unpaired_file2}" are the output files.

# Lastly, FastQC is run on the trimmed paired files:
apptainer exec /containers/apptainer/fastqc-0.12.1.sif fastqc -o ${output_dir} -t 4 ${trimmed_file1} ${trimmed_file2}
# CPU threads per task is 4 (-t 4), same as trimmomatic.