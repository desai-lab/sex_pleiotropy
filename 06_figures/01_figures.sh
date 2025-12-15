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

python fig2_figS4.py

python figure3_processing.py
python figure4_processing.py

module load python
source activate gge_R_figures
module unload python

Rscript fig1.R
Rscript figS2.R
Rscript figS3.R
Rscript fig3_figS5.R
Rscript fig4_figS6.R