## Correcting the original reference genome by pooling data from 65 clones 

6 clones that were removed from all further analyses due to suspected whole-genome duplications, leaving a total of 65 clones.  

Combine these alignments and polish the original genome with Pilon to correct for any ancestral mutations:
```  
sbatch 01_gge_ref_corr_pipeline_65_clones.sh
```  

Calculate coverage across genome for all 65 clones:
```
sbatch 02_calculate_coverage.sh
```