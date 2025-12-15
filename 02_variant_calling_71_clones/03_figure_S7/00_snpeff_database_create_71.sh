#Here's the code I ran in an interactive session

# Load necessary modules and activate environment
module load python
source activate gge_figures

cd /n/home03/spai/.conda/envs/gge_figures/share/snpeff-5.2-1

mkdir data
cd data

mkdir w303_gge_corr_71
cd w303_gge_corr_71

cp /n/holystore01/LABS/desai_lab/Lab/spai/GGE/GGE_pipeline_final/01_reference_correction_71_clones/w303_gge_corr_71.fasta sequences.fa
cp /n/holystore01/LABS/desai_lab/Lab/spai/GGE/GGE_pipeline_final/01_reference_correction_71_clones/w303_gge_corr_71.gff genes.gff

cd ../..

nano snpEff.config
# added this line to the end of the .config (by pressing Esc + / in nano):
# # Saccharomyces cerevisiae w303 corrected 71
# w303_gge_corr_71.genome : w303_gge_corr_71

java -jar snpEff.jar build -gff3 -noCheckCds -noCheckProtein -v w303_gge_corr_71