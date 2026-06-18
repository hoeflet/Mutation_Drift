# Mutation_Drift
Data analysis and bioinformatics pipelines used for the project "Mutational Input, Genetic Drift and Defective Viral Genomes Shape Genetic Diversity in Herpes Simplex Virus Type 1".

## NGS Pipelines
NGS Analysis was performed on the Freie Universität Berlin HPC cluster Curta. Raw data (.fastq) was trimmed with Trimmomatic, mapped to the reference with BWA and variant called with Freebayes, LoFreq and BCFtools. The full pipeline can be found here (```NGS_pipeline.py```). Recombination mapping was performed with ViReMa using ```NGS_recomb.py```.

## Data Analysis
All scripts used to analyze and plot data is given within the ```MutDrift.ipynb``` jupyter notebook. Cells are annotated with markdown cells above.

## Reference

Please cite our paper if you reuse any of our analysis or pipelines.

*Paper in preparation, citation will be updated*
