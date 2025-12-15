## Variant Calling of 71 clones with Corrected Reference

Perform variant calling for all 71 clones:
```  
sbatch 01_gge_pipeline_corr_71_clones.sh
sbatch 02_gge_pipeline_corr_database_71_clones.sh
```

In an interactive session, create a SnpEff database following steps in `00_snpeff_database_create_71.sh`.  

Process variants and make figure S7 - this shows the 6 clones that were removed from all further analyses due to suspected whole-genome duplications:
```
cd 03_figure_S7
sbatch 01_gge_vcf_processing_71_clones.sh
```