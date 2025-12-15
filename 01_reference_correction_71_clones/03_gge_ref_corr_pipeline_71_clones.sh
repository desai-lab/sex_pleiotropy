#!/bin/bash
#SBATCH -c 1                              # Number of cores
#SBATCH -t 0-12:00                        # Runtime in D-HH:MM
#SBATCH -p serial_requeue                 # Partition to submit to
#SBATCH --mem-per-cpu=240000                       # Memory pool for all cores
#SBATCH -o shell_outputs/myoutput_%j.out      # File to which STDOUT will be written
#SBATCH -e shell_outputs/myerrors_%j.err      # File to which STDERR will be written
#SBATCH --open-mode=append
#SBATCH --mail-type=ALL                   # Email notification type
#SBATCH --mail-user=shreyaspai@g.harvard.edu # Email for notifications

module load python
source activate gge_variant_calling
module unload python

mkdir -p shell_outputs logs

# Merge BAM files
samtools merge combined_ref.bam ./pre_corr_outputs/final_bam/*final.bam

# Sort and index the merged BAM file using Picard
picard SortSam \
    INPUT=combined_ref.bam \
    OUTPUT=combined_ref.final.bam \
    SORT_ORDER=coordinate \
    CREATE_INDEX=true \
    2> logs/combined_ref_sorting.log

# Validate the sorted BAM file
picard ValidateSamFile \
    INPUT=combined_ref.final.bam \
    OUTPUT=combined_ref.validate.txt \
    MODE=SUMMARY \
    2> logs/combined_ref_validate.log

# Run Pilon for genome refinement
pilon \
    --genome ../uncorrected_reference/w303_ref.fasta \
    --bam combined_ref.final.bam \
    --output w303_gge_corr_71 \
    --outdir pilon_output/ \
    --changes \
    --tracks \
    -Xmx240g

# Update GFF and fix FASTA headers
python update_gff_and_fix_fasta_headers_71.py

# Index the reference genome with samtools
samtools faidx w303_gge_corr_71.fasta

# Index the reference genome with BWA
# This will create 5 additional files: .amb, .ann, .bwt, .pac, .sa
bwa index -p w303_gge_corr_71 w303_gge_corr_71.fasta

# Create a sequence dictionary using Picard
picard CreateSequenceDictionary \
    R=w303_gge_corr_71.fasta \
    O=w303_gge_corr_71.dict
