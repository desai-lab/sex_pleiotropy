#!/bin/bash
#SBATCH -c 6                              # Number of cores
#SBATCH -t 0-03:00                        # Runtime in D-HH:MM
#SBATCH -p serial_requeue                 # Partition to submit to
#SBATCH --mem=64000                       # Memory pool for all cores
#SBATCH -o shell_outputs/myoutput_%j.out      # File to which STDOUT will be written
#SBATCH -e shell_outputs/myerrors_%j.err      # File to which STDERR will be written
#SBATCH --open-mode=append
#SBATCH --mail-type=ALL                   # Email notification type
#SBATCH --mail-user=shreyaspai@g.harvard.edu # Email for notifications

module load python
source activate gge_arc

reference="../01_reference_correction_71_clones/w303_gge_corr_71.fasta"

# 8B) Generating a database from the g.vcf files
# Prepare the list of -V arguments for GenomicsDBImport
vcfs=$(find *.g.vcf | awk '{print "-V " $0}' | tr '\n' ' ')

gatk GenomicsDBImport \
    -R "${reference}" \
    --genomicsdb-workspace-path database_gge \
    -L chrI -L chrII -L chrIII -L chrIV -L chrV \
    -L chrVI -L chrVII -L chrVIII -L chrIX -L chrX \
    -L chrXI -L chrXII -L chrXIII -L chrXIV -L chrXV \
    -L chrXVI -L chrMito -L 2-micron \
    ${vcfs} \
    --reader-threads 6

# 8C) Genotyping from the database
gatk GenotypeGVCFs \
    -R ${reference} \
    -V gendb://database_gge \
    -O GGE_clones_corr.vcf \
    --heterozygosity 0.005 \
    2> genotyping.log \
    --sample-ploidy 1

# Consolidation Step: Organize all outputs
output_dir="./organized_outputs"
mkdir -p ${output_dir}/{trimmed_fastq,sam_files,sorted_bam,dedup_bam,final_bam,metrics,gvcf_files,logs,vcf}
find . -name "*.trimmed_*" -exec mv {} ${output_dir}/trimmed_fastq/ \;
find . -name "*.sam" -exec mv {} ${output_dir}/sam_files/ \;
find . -name "*.sorted.*" -exec mv {} ${output_dir}/sorted_bam/ \;
find . -name "*.dedup.bam" -exec mv {} ${output_dir}/dedup_bam/ \;
find . -name "*.dedup_metrics.txt" -exec mv {} ${output_dir}/metrics/ \;
find . -name "*.final.*" -exec mv {} ${output_dir}/final_bam/ \;
find . -name "*.validate.txt" -exec mv {} ${output_dir}/metrics/ \;
find . -name "*.g.vcf*" -exec mv {} ${output_dir}/gvcf_files/ \;
find . -name "*.log" -exec mv {} ${output_dir}/logs/ \;
mv GGE_clones_corr.vcf* ${output_dir}/vcf/
mv genotyping.log ${output_dir}/logs/
mv database_gge ${output_dir}/
find shell_outputs -type f -name "*.out" -exec mv {} ${output_dir}/logs/ \;
find shell_outputs -type f -name "*.err" -exec mv {} ${output_dir}/logs/ \;