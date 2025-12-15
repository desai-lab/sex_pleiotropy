#!/bin/bash
#SBATCH -c 1                              # Number of cores
#SBATCH -t 0-04:00                        # Runtime in D-HH:MM
#SBATCH -p serial_requeue                 # Partition to submit to
#SBATCH --mem-per-cpu=20000                      # Memory
#SBATCH -o shell_outputs/myoutput_%j.out      # File to which STDOUT will be written
#SBATCH -e shell_outputs/myerrors_%j.err      # File to which STDERR will be written
#SBATCH --open-mode=append
#SBATCH --mail-type=ALL                   # Email notification type
#SBATCH --mail-user=shreyaspai@g.harvard.edu # Email for notifications

module load python
source activate gge_variant_calling

bedtools genomecov -d -ibam combined_ref.final.bam > coverage.txt
