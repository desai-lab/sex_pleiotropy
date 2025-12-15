## Correcting the original reference genome by pooling data from all 71 clones 

For all 71 clones, align reads to original reference genome (in parallel):
```  
sbatch 01_gge_pre_corr_pipeline.sh
```

Consolidate the files:
```  
sbatch 02_consolidate.sh
```  

Combine these alignments and polish the original genome with Pilon to correct for any ancestral mutations:
```  
sbatch 03_gge_ref_corr_pipeline_71_clones.sh
```  

Calculate coverage across genome for all 71 clones:
```
sbatch 04_calculate_coverage.sh
```