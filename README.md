# Pacific oyster panel design and testing
Code to support the design and testing an amplicon panel for Pacific oyster (_Crassostrea gigas_), designed for the use of the authors and associated [manuscript](#toadd), and comes with no guarantees for usefulness for any other purposes.       

Requires Linux or Mac OS, and all shell scripts are run from the main directory.    

Skip to [panel design section](https://github.com/bensutherland/ms_oyster_panel/blob/main/README.md#panel-design)       
Skip to [chromosome coordinates section](https://github.com/bensutherland/ms_oyster_panel/blob/main/README.md#chromosome-positions)       
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


## Additional Analyses for Pilot Study
Use the script `01_scripts/use_multiple_run_data.R` to analyze data from the two different runs that were sequenced as part of the pilot study. This script will allow you to choose the best technical replicate (i.e., the one with the most typed markers) and combine into a single genind file for analysis. It will also evaluate technical replicates within a sequencing run.        
Once the best individual per technical reps is identified, move to `01_scripts/simple_pop_stats_Cgig_analysis_*.R`, which will allow you to analyze the data (requires the [simple_pop_stats](https://github.com/bensutherland/simple_pop_stats) repo).     



## Chromosome positions ##


## Panel testing ##

