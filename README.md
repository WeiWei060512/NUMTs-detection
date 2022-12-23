# NUMTs-detection
This repository contains scripts for data processing and analysis used in the publication:

Wei W, Schon K, Elgar G, Orioli A, Tanguy M, Giess A, Tischkowitz M, Caulfield M, Chinnery PF. Nuclear-embedded mitochondrial DNA sequences in 66,083 human genomes. Nature 611, 105–114 (2022). https://doi.org/10.1038/s41586-022-05288-7 

Note: It's not a software to do general analysis and won't be maintained. If you do use the code in your publications, please cite our paper.

1. NUMTs and breakpoints detection 

NUMTs_detection.sh

searchBreakpoint_fromblatoutputs.py

searchNumtCluster_fromDiscordantReads.py

groupNumtCluster_fromMultipleSamples.py


2. enrichment analysis

enrichment_creatingRefgenome.py

enrichment_simulation.py

3. mtDNA variants calling

mtVariantCalling.sh

mtVariantCalling_MToolBox.conf


4. NUMTs methylation detection

nanopolish_methylationDetection.sh

5. NUMT variants calling

VarDetection_fromDiscSplitReads.sh

generateVariantTable.Human.py

generateVariantTable.HumanChimp.py

6. Circos plots

circos_allNUMTs.conf

confs

## Data Availability
Whole genome sequence data that support the findings of this study can be analysed on the Genomics England data warehouse through https://www.genomicsengland.co.uk/understanding-genomics/data/

