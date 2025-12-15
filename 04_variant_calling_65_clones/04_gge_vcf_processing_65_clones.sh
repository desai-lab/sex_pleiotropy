#!/bin/bash
#SBATCH -c 1                              # Number of cores
#SBATCH -t 0-00:20                        # Runtime in D-HH:MM
#SBATCH -p serial_requeue                 # Partition to submit to
#SBATCH --mem-per-cpu=20000                # Memory per core
#SBATCH -o shell_outputs/myoutput_%A_%a.out # File to which STDOUT will be written
#SBATCH -e shell_outputs/myerrors_%A_%a.err # File to which STDERR will be written
#SBATCH --open-mode=append
#SBATCH --mail-type=ALL                   # Email notification type
#SBATCH --mail-user=shreyaspai@g.harvard.edu # Email for notifications

# Load necessary modules and activate environment
module load python
source activate gge_figures
module unload python

# Define file paths
CLONES_VCF="./organized_outputs/vcf/GGE_clones_corr.vcf"

# See snpeff_database_create_65.sh for how I created the appropriate SnpEff database for w303_gge_corr_65

# Annotating variants with snpEff
snpEff w303_gge_corr_65 $CLONES_VCF > GGE_clones_corr_ann.vcf

# Converting VCF to TSV for further processing
gatk VariantsToTable \
    -V GGE_clones_corr_ann.vcf \
    -F CHROM -F POS -F REF -F ALT -F QUAL -F ANN -GF AD \
    -O GGE_clones_corr_ann.tsv

# Split multi-allelic records containing * into individual records
# Do not -- Delete all other multi-allelic records (likely to be erroneous positions)
# also filter out alleles where ALT = * (spanning/overlapping deletion)
python split_multi_allelic_records.py GGE_clones_corr_ann.tsv

# Filter out telomeric variants
# Filter out mitochondrial variants
# Filter out variants that appear in >=2 populations (likely to be in ancestor)
# Filter out variants from extremely low coverage regions (<10x average coverage per clone)
# Filter out variants with any clone's alt allele frequency < 0.75
# Filter out variants with low alt allele counts (in all clones they appear) & set alt counts <5 to 0
python filter_variants.py GGE_clones_corr_ann_split.tsv

# Then simplify annotations (VLTE-style)
python simplify_annotations.py GGE_clones_corr_ann_split_filtered.tsv

# Group nearby variants (<25bp) into "mutation groups"
python group_variants.py GGE_clones_corr_ann_split_filtered_simplified.tsv

# Substitute clone_ids for filenames which currently serve as headers for each clone
python change_headers.py GGE_clones_corr_ann_split_filtered_simplified_grouped.tsv