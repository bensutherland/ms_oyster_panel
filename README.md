# ms_oyster_panel
Repository for code to support the manuscript to design a Pacific oyster amplicon panel

Note: this repository is for the designed use of the author only and comes with no guarantees of usefulness for anything else.       

Data inputs:     
- Single SNP data in genepop format from Sutherland et al. 2020 (Evol. Appl.)        
- Reference genome for Pacific oyster used in identifying markers (Zhang et al. 2012; Sutherland et al. 2020)        
- Single SNP data populations output in VCF format

Note: the genepop has been filtered for -r 0.7 (70% of individuals) in all of 15 populations, with a --min-maf >= 0.01 with HWE outliers removed (as per Sutherland et al. 2020).        



