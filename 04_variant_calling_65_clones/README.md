## Variant Calling of 65 clones with Corrected Reference

Perform variant calling for 65 clones with newly corrected reference genome:
```  
sbatch 01_gge_pipeline_corr_65_clones.sh
sbatch 02_gge_pipeline_corr_database_65_clones.sh
```

In an interactive session, create a SnpEff database following steps in `03_snpeff_database_create_65.sh`.  

Process variants to get final table `GGE_clones_final.tsv`:
```
sbatch 04_gge_vcf_processing_65_clones.sh
```