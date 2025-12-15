#!/bin/bash
#SBATCH -c 1                              # Number of cores
#SBATCH -t 0-02:00                        # Runtime in D-HH:MM
#SBATCH -p serial_requeue                 # Partition to submit to
#SBATCH --mem-per-cpu=20000                # Memory per core
#SBATCH -o shell_outputs/myoutput_%A_%a.out # File to which STDOUT will be written
#SBATCH -e shell_outputs/myerrors_%A_%a.err # File to which STDERR will be written
#SBATCH --array=0-70                      # Job array range (based on sample count)
#SBATCH --open-mode=append
#SBATCH --mail-type=ALL                   # Email notification type
#SBATCH --mail-user=shreyaspai@g.harvard.edu # Email for notifications

module load python
source activate gge_variant_calling

mkdir -p shell_outputs

# Define the path to the sample list
sample_list_file="../gge_clones_fastq_files_list.txt"

# Extract the R1 file based on the job array index
r1_file=$(sed -n "$((2 * SLURM_ARRAY_TASK_ID + 1))p" ${sample_list_file})
r2_file=$(sed -n "$((2 * SLURM_ARRAY_TASK_ID + 2))p" ${sample_list_file})

# Get the base name for output files by removing .R1/.R2.fastq.gz suffix
sample_base=$(basename ${r1_file} .R1.fastq.gz)

# Define the directory for FASTQ files and reference genome path
fastq_dir="../../gge_clones_fastq_files"
reference="../uncorrected_reference/w303_ref.fasta"

# Step 1: Read Trimming
NGmerge -1 ${fastq_dir}/${sample_base}.R1.fastq.gz \
        -2 ${fastq_dir}/${sample_base}.R2.fastq.gz \
        -a -v -o ${sample_base}.trimmed

# Step 2: Alignment Against Reference Genome

# #Indexing the reference genome for bwa -- do this only once
# bwa index -p w303_ref ${reference}

bwa mem -M -t 1 -R "@RG\tID:HTLMC.1\tSM:${sample_base}\tPL:ILLUMINA" ../uncorrected_reference/w303_ref \
        ${sample_base}.trimmed_1.fastq.gz ${sample_base}.trimmed_2.fastq.gz > ${sample_base}.sam \
        2> ${sample_base}.bwa.log

# Step 3: Conversion from SAM to BAM
picard SortSam \
    I=${sample_base}.sam \
    O=${sample_base}.sorted.bam \
    SORT_ORDER=coordinate \
    CREATE_INDEX=true

# Step 4: Marking Duplicates
picard MarkDuplicates \
    I=${sample_base}.sorted.bam \
    O=${sample_base}.dedup.bam \
    METRICS_FILE=${sample_base}.dedup_metrics.txt \
    REMOVE_DUPLICATES=false \
    TAGGING_POLICY=All 2>> dedup.log

# Step 5: Resorting and Reindexing
picard SortSam \
    I=${sample_base}.dedup.bam \
    O=${sample_base}.final.bam \
    SORT_ORDER=coordinate \
    CREATE_INDEX=true 2>> final_sorting.log

# Step 6: Validating BAM Files
picard ValidateSamFile \
    I=${sample_base}.final.bam \
    O=${sample_base}.validate.txt \
    MODE=SUMMARY
