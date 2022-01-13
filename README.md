# ms_oyster_panel
Repository for code to support the manuscript to design a Pacific oyster amplicon panel

Note: this repository is for the designed use of the author only and comes with no guarantees of usefulness for anything else.       

Data inputs:     
- Single SNP data in plink format from Sutherland et al. 2020 (Evol. Appl.)        
- Single SNP data populations output in VCF format       
- Reference genome for Pacific oyster used in identifying markers (Zhang et al. 2012; Sutherland et al. 2020)         
Download this version, GenBank format: https://www.ncbi.nlm.nih.gov/assembly/GCF_000297895.1/        

Note: the genepop has been filtered for -r 0.7 (70% of individuals) in all of 15 populations, with a --min-maf >= 0.01 with HWE outliers removed (as per Sutherland et al. 2020).        

### 00. Filter and characterize variants ###
Use the script `01_scripts/01_identifying_markers.R` run interactively in R. This script will do the following:      
1. Load packages and set working directory (need to change if on non-design system)
2. Import data and select only the BC naturalized samples
3. Perform minor allele frequency (MAF) filter to remove low MAF variants
4. Calculate per locus stats observed heterozygosity (Hobs) and FST (averaged per locus)
5. Generate comparative plots
6. Output a csv of marker per locus stats for all MAF-filtered variants as `03_marker_selection/all_markers_per_locus_stats_incl_hobs_MAF.csv`     

