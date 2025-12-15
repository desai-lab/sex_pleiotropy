#!/bin/bash
#SBATCH -c 1                              # Number of cores
#SBATCH -t 0-00:01                        # Runtime in D-HH:MM
#SBATCH -p serial_requeue                 # Partition to submit to
#SBATCH --mem-per-cpu=1000                # Memory per core
#SBATCH -o shell_outputs/myoutput_%A_%a.out # File to which STDOUT will be written
#SBATCH -e shell_outputs/myerrors_%A_%a.err # File to which STDERR will be written
#SBATCH --open-mode=append
#SBATCH --mail-type=ALL                   # Email notification type
#SBATCH --mail-user=shreyaspai@g.harvard.edu # Email for notifications

# Consolidation Step: Organize all outputs
output_dir="./pre_corr_outputs"
mkdir -p ${output_dir}/{trimmed_fastq,sam_files,sorted_bam,dedup_bam,final_bam,metrics,gvcf_files,logs,vcf}
find . -name "*.trimmed_*" -exec mv {} ${output_dir}/trimmed_fastq/ \;
find . -name "*.sam" -exec mv {} ${output_dir}/sam_files/ \;
find . -name "*.sorted.*" -exec mv {} ${output_dir}/sorted_bam/ \;
find . -name "*.dedup.bam" -exec mv {} ${output_dir}/dedup_bam/ \;
find . -name "*.dedup_metrics.txt" -exec mv {} ${output_dir}/metrics/ \;
find . -name "*.final.*" -exec mv {} ${output_dir}/final_bam/ \;
find . -name "*.validate.txt" -exec mv {} ${output_dir}/metrics/ \;
find . -name "*.log" -exec mv {} ${output_dir}/logs/ \;
find shell_outputs -type f -name "*.out" -exec mv {} ${output_dir}/logs/ \;
find shell_outputs -type f -name "*.err" -exec mv {} ${output_dir}/logs/ \;