# Pacific oyster panel design and testing
Code to support the design and testing an amplicon panel for Pacific oyster (_Crassostrea gigas_), designed for the use of the authors and associated [manuscript](https://www.biorxiv.org/content/10.1101/2023.08.26.554329v1), and comes with no guarantees for usefulness for any other purposes.       

Requires Linux or Mac OS, and all shell scripts are run from the main directory.    

Skip to [panel design section](https://github.com/bensutherland/ms_oyster_panel/blob/main/README.md#panel-design)       
Skip to [chromosome positions section](https://github.com/bensutherland/ms_oyster_panel/blob/main/README.md#chromosome-positions)       
Skip to [panel testing section](https://github.com/bensutherland/ms_oyster_panel/blob/main/README.md#panel-testing)         


## Panel design ##
#### Requirements:      
- bedtools     
- R 
- Eric Normandeau's [scripts](https://github.com/enormandeau/Scripts)        
 

#### Data inputs:     
- Prefiltered single SNP data in plink format from Sutherland et al. 2020 (Evol. Appl.) (here) #TODO       
- Prefiltered single SNP data populations output in VCF format (here) #TODO      
- Reference genome for Pacific oyster used in identifying markers (Zhang et al. 2012; Sutherland et al. 2020) [download GenBank version here](https://www.ncbi.nlm.nih.gov/assembly/GCF_000297895.1/)         


### 01. Filter and characterize variants ###
Interactively run `01_scripts/01_identifying_markers.R` to       
1. Import data and select only the British Columbia (BC) naturalized samples
2. Perform minor allele frequency (MAF) filter to remove variants under MAF in BC populations
3. Calculate per locus stats observed heterozygosity (Hobs) and FST (averaged per locus)
4. Generate comparative per locus stat plots

**Outputs**       
- `03_marker_selection/all_markers_per_locus_stats_incl_hobs_MAF.csv`     
e.g.,     
```
mname,Fit,Fst,Fis,maf.vec,Hobs
100026_37_C,-0.011,0.00777,-0.0192,0.0173,0.0347 
```


### 02. Select custom markers (private alleles)
Custom markers will be included in the panel, specifically private allele variants specific to cultured BC populations. This will use code and resources adapted from `https://github.com/bensutherland/ms_oyster_popgen`, specifically `00_archive/my_cols.csv` and `01_scripts/private_alleles.r`.      

Interactively run `01_scripts/private_alleles.r` to identify marker names (mnames) of high frequency private alleles from DPB and GUR populations.    

**Outputs:**    
- `03_marker_selection/per_repunit_private_allele_tally_all_data.csv`     
- `03_marker_selection/DPB_selected_PA_mnames.csv`    
- `03_marker_selection/GUR_selected_PA_mnames.csv`    


### 03. Select markers (high FST or high HOBS)
Interactively run `01_scripts/01b_collect_mnames.sh` and set the number of markers that you want to retain based on high HOBS and high FST characteristics.     

**Outputs**
- `04_extract_loci/selected_mnames.csv`        
e.g.,     
```
504898_23_C, top_FST
375205_31_T, top_FST
```


### 04. Obtain info from VCF on selected markers ###
Interactively run `01_scripts/02_info_from_vcf.sh` to obtain specific information about the markers from the larger VCF       

**Outputs** 
- `04_extract_loci/vcf_selection.csv`             
This is a text file with fields (1) chr; (2) pos of SNP in ref genome; (3) info about marker; (4) ref allele; (5) alt allele, e.g.,                 
```
JH816256.1,99460,100388:26:-,G,A
```

### 05. Prepare BED file for extracting selected marker sequence
Interactively run `01_scripts/03_prepare_bed_file.sh` to prepare a bed file based on the selected markers

**Outputs**
- `04_extract_loci/vcf_selection.bed`
e.g., 
```
#TODO note: this will include a fourth column that can be used for matching the output fasta back to the VCF
```


### 06. Extract FASTA from reference genome
Run `01_scripts/04_extract_from_reference.sh` to extract selected marker flanking sequence from the genome

**Outputs**     
- `04_extract_loci/vcf_selection.fa`        
- `04_extract_loci/selected_chr_and_seq.txt` (tab delim version)    


### 07. Create design submission file
Interactively run `01_scripts/05_make_submission_form.R` to connect sequence data to selected marker information from the VCF.         

**Outputs**        
- `05_submission_form/seq_and_minfo_all_data.csv` (full information for data checking)
- `05_submission_form/seq_and_minfo_for_submission.csv` (submission info only)
The submission csv has fields (1) marker name; (2) chr; (3) ref allele (based on genome); (4) alt allele; (5) strand; (6) marker type; (7) priority level; (8) formatted seq (e.g., ATGC[A/G]ATGC). More details are available in the script.      

#TODO# fix: this section 
Note that the reference allele is given that designation based on the reference genome nucleotide and doesn't always indicate in that manner in the VCF. A switch will occur in the Rscript to switch these identities, if the VCF alt allele is the allele in the reference genome. 
#END TODO# fix: this section


### 08. Additional data checking (optional)
Interactively run `01_scripts/confirm_FST.R` to confirm that the selected markers for high FST are indeed leading to higher population-level FST estimates.        

**Inputs**
- `03_marker_selection/adegenet_output.RData` (generated above)
- `04_extract_loci/selected_mnames.csv` 


#### Also see information [here](#TODO) for additional information on the bed and VCF formats for extraction.       


[simple_pop_stats](https://github.com/bensutherland/simple_pop_stats) repo).     



## Chromosome positions ##
Currently implemented in an interactive Rscript `100-bowtie2-amplicon-mapping-roslin-genome.Rmd` [here](https://github.com/bensutherland/ms_oyster_panel/tree/main/100_amplicon_mapping) 

#### Requirements
- bowtie2
- samtools
- R

#### Data inputs
- Marker information file `thermo_submitted_sutherland_cgig_v.0.1_2022-02-08_marker_data_only_2023-02-21.xlsx` available as Additional File X [here](#toadd)      
- Reference genome (Roslin Institute) available [here](https://www.ncbi.nlm.nih.gov/assembly/GCF_902806645.1)      

**Outputs**
- Alignment sam file
- Plot of distribution of markers across the reference genome chromosomes including singly or multiply aligning markers


## Panel testing ##
#### Requirements:      
- [amplitools](https://github.com/bensutherland/amplitools)       
- [simple_pop_stats](https://github.com/bensutherland/simple_pop_stats)      
- [R](https://www.r-project.org/)             

Create a folder that contains both `amplitools` and `simple_pop_stats`, as both of these repositories will be used to analyze the data.       


#### Data inputs:     
- VariantCaller output files (tab-delim but .xls suffix) from [FigShare](https://doi.org/10.6084/m9.figshare.23646471.v1)       
Specifically, the manuscript analysis will use `R_2022_08_04_S5XL.xls` and `R_2022_10_07_S5XL.xls`.       

### 01. Load results from available amplicon sequencing VariantCaller outputs 
Copy VariantCaller output files into the `amplitools/02_input_data` folder.     

Open `amplitools/01_scripts/00_initiator.R` in Rstudio and source the script. This will load amplitools.       

Open `amplitools/01_scripts/demo_analysis.R` for a demonstration analysis (also shown in brief below).      

From within R, convert proton results to genepop results using the following amplitools script:      
```
proton_to_genepop(hotspot_only=TRUE, neg_control="BLANK")         
# hotspot_only (T/F) indicates whether only hotspot variants will be considered, or all variants including novel variants
# neg_control is a string found in only negative control samples to filter them out
```

Then from the amplitools main directory, finalize the genepop format using the following for each file:      
`amplitools/01_scripts/format_genepop.sh <filename>`        

**Outputs**
- `*.gen` files prepared in `amplitools/02_input_data/prepped_matrices`     


### 02. Compare technical replicate samples and retain the best replicates 
Copy the output `*.gen` files from above into `simple_pop_stats/02_input_data`.        

In Rstudio clear the workspace, then source `simple_pop_stats/01_scripts/simple_pop_stats_start.R`     
Also source the development script `simple_pop_stats/01_scripts/dev/comp_tech_reps.R`       

Run the following:      
`comp_tech_reps(format_type = "amplitools", max_missing = 0.5)`        
Where `max_missing` indicates the maximum missing data to retain a sample for the technical replicate comparison.       

Save out the produced genind object, which contains only the best of the technical replicate samples:      
`save(obj_nr_best, file = "02_input_data/obj_nr_best_<date>.RData")`      

Note: if you need to restart at any future time, you can always reload this file by:     
`load("02_input_data/obj_nr_best_<date>.RData")`      

**Outputs**      
- obj_nr_best , a single sample per tech rep for population genetic analysis (below)
- #TODO add other items


### 03. General population genetic characterization of pilot study
With `simple_pop_stats` functions still sourced (i.e., not cleared space), in Rstudio, run interactively `ms_oyster_panel/01_scripts/sps_popgen_analysis.R`        
note: if you need to reload the `obj_nr_best` see above       
note: if you are not running tech replicates, can use `load_genepop(datatype="SNP")` to select your genepop

Follow the instructions of the script to:     
1. Prepare data including adding population attribute to the genepop
2. Characterize missing data and filter as needed
3. Remove monomorphic loci
4. Generate per-locus statistics
5. Run multivariate and differentiation analyses
6. Identify private alleles at the regional or population level
7. Convert the genepop to a rubias file for downstream use in parentage (see below)

**Outputs**
- various simple_pop_stats output in 03_results
- a rubias file will be output to `amplitools/03_prepped_data/cgig_all_rubias.txt`    

### 04. Parentage analysis
Requires that the rubias file was produced in the preceding step. This will contain three generations of samples, `VIU_F0`, `VIU_F1`, `VIU_F2` for parentage analysis.      
Clear the workspace, then source `amplitools/00_initiator.R`      

Estimate log likelihoods from the data and simulated siblings and parents, then calculate on your existing data:     
Run the following two scripts to analyze each set of contrasts:       
```
# For F1 vs. F2 (parents vs. offspring)
ckmr_from_rubias(input.FN = "03_prepped_data/cgig_all_rubias.txt", parent_pop = "VIU_F1", offspring_pop = "VIU_F2", cutoff = 5)

# For F0 vs. F1 (grandparents vs. parents)
ckmr_from_rubias(input.FN = "03_prepped_data/cgig_all_rubias.txt", parent_pop = "VIU_F0", offspring_pop = "VIU_F1", cutoff = 5)

```

**Outputs**
- Figure of logl from simulated relationships
- Parent-offspring relationships above the cutoff
- Parent and offspring full-sib results above the cutoff
- Retained individuals for both parent and offspring populations



