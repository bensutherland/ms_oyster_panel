# ms_oyster_panel
Code to support the design of a Pacific oyster amplicon panel. Starts with a VCF and plink data export from Stacks populations module, and uses reference genome, to produce a submission-ready text file for amplicon design.      

_Note: this repository is for the designed use of the author only and comes with no guarantees of usefulness for anything else._       

Requires basic Linux or Mac OS.      
All shell scripts are run from the main directory.     

#### Requirements:      
- Input files (see below)
- bedtools     
- [Eric Normandeau's scripts](https://github.com/enormandeau/Scripts)        
- R packages as designated in the Rscripts     
 

#### Data inputs:     
- Single SNP data in plink format from Sutherland et al. 2020 (Evol. Appl.)        
- Single SNP data populations output in VCF format       
- Reference genome for Pacific oyster used in identifying markers (Zhang et al. 2012; Sutherland et al. 2020)         
Download this version, GenBank format: https://www.ncbi.nlm.nih.gov/assembly/GCF_000297895.1/        

Note: in Sutherland et al. 2020, the genetic data was been filtered to retain markers present in 70% of individuals in all of 15 populations, with a --min-maf >= 0.01, and with HWE outliers removed.        


### 01. Filter and characterize variants ###
Use the script `01_scripts/01_identifying_markers.R` run interactively in R. This script will do the following:      
1. Import data and select only the BC naturalized samples
2. Perform minor allele frequency (MAF) filter to remove low MAF variants wrt. BC pops
3. Calculate per locus stats observed heterozygosity (Hobs) and FST (averaged per locus)
4. Generate comparative per locus stats plots

...output per locus stats for filtered variants:      
`03_marker_selection/all_markers_per_locus_stats_incl_hobs_MAF.csv`     

...with the following format:    
```
mname,Fit,Fst,Fis,maf.vec,Hobs
100026_37_C,-0.011,0.00777,-0.0192,0.0173,0.0347 
```


### 02. Additional markers to include
Custom markers will be included in the panel. Currently this includes markers that are private alleles for a couple of populations in culture only in British Columbia. This will use code and resources adapted from `https://github.com/bensutherland/ms_oyster_popgen`, specifically the following:     
`00_archive/my_cols.csv`         
`01_scripts/private_alleles.r`      
...which have been named with the same filenames/ folders in the current repo, but updated mainly to allow output of specific marker names.     

Run `01_scripts/private_alleles.r`      
...this will identify mnames of high frequency private alleles from DPB and GUR populations.    
Outputs:    
`03_marker_selection/per_repunit_private_allele_tally_all_data.csv`     
`03_marker_selection/DPB_selected_PA_mnames.csv`    
`03_marker_selection/GUR_selected_PA_mnames.csv`    


### 03. Select the top markers to create whitelists ###
Set the number of markers you want to collect from each, then run `01_scripts/01b_collect_mnames.sh`    
Output will be rows of mnames and mtype, no header.      

`04_extract_loci/selected_mnames.csv`        
```
504898_23_C, top_FST
375205_31_T, top_FST
563276_7_A, top_FST
355586_10_A, top_FST
```


### 04. Collect the whitelist marker info from VCF ###
Obtain specific information about these markers from the larger VCF:       
`01_scripts/02_info_from_vcf.sh`          
...outputs `04_extract_loci/vcf_selection.csv`, which contains (for example):             


```
JH816256.1,99460,100388:26:-,G,A

# where the fields are: 
# (1) chromosome
# (2) position of the SNP in the reference genome
# (3) information about the marker
# (4) reference allele
# (5) variant allele

```


### 05. Prepare BED file
Use the following to prepare a bed file from the relevant lines of the vcf
`01_scripts/03_prepare_bed_file.sh`     
...will produce `04_extract_loci/vcf_selection.bed`
note: this will include a fourth column that can be used for matching the output fasta back to the VCF (i.e., mname).     



### 06. Extract FASTA from reference genome
Then use the following to extract the relevant sequence from the genome
`01_scripts/04_extract_from_reference.sh`       
...will produce `04_extract_loci/vcf_selection.fa`        


### 7. Bring all back together and create submission file
Use the following Rscript interactively to join the sequence data (produced from the fasta) to the marker information (produced from the VCF):        
`05_make_submission_form.R`       

This script will do the following:     

...this will produce `05_submission_form/seq_and_minfo_for_submission.csv`, which contains the following:    
1. marker name;     
2. chromosome (scaffold);     
3. reference allele (reference based on the reference genome)   
4. alternate allele
5. strand;     
6. marker type;    
7. priority level;     
8. formatted sequence (...ATGC[A/G]ATGC...)     

Note that the reference allele is given that designation based on the reference genome nucleotide and doesn't always indicate in that manner in the VCF. A switch will occur in the Rscript to switch these identities, if the VCF alt allele is the allele in the reference genome. 

For data checking purposes, a full dataframe as generated through the Rscript process is also output as `05_submission_form/seq_and_minfo_all_data.csv`       


#### Now the data can be submitted for marker design ####


### Additional information on VCF formats ####
This section is for review only, but contains some relevant information about formats.     
To extract from the reference genome, we will use bedtools combined with a bedfile.     
The bedfile will be in the format of:     
```
chr1	5	10	<VCF_info_field>
```
...a tab-delimited file with positions refering to select positions in the contig/scaffold starting counting from position 0.      

Get a better view of an accession in the fasta. First, use `enormandeau/Scripts/fasta_unrwap.py` to get a better view of the fasta file. 
`fasta_unwrap.py 02_input_data/GCA_000297895.1_oyster_v9_genomic.fna 02_input_data/GCA_000297895.1_oyster_v9_genomic_unwrapped.fna`     

Take an example from the VCF information:       
`JH816256.1,99460,100388:26:-,G,A`       

Review the fasta consensus sequence for this from the stacks fasta file:      
```
>CLocus_100388 [JH816256.1, 99396, -]
TGCATCTGTCCTCTTTTCTGTGCTTCATTGATGGTCATTTTGGATACCTGTTTGACTTAACTGATTTAGATAAGATGATCATGTGTTGTG
```

Extract from the reference genome:       
`grep -A1 JH816256.1 02_input_data/GCA_000297895.1_oyster_v9_genomic_unwrapped.fna > 02_input_data/example_record_JH816256.1.txt`     

Looking at the reference genome record (open in text editor) you will see that the position 99396 (referred to in the stacks fasta accession) indicates the starting character of the reverse complement of the stacks fasta record.     
Further, looking at position 99460 in the ref genome, you will see that this is the 'G' nucleotide, which is the reference allele (i.e., TGAAGC). Note: this would be in a 1-based count.        

We are therefore looking to extend 200 bp on either side of this variant (variant being at 99459-99460 in 0th positional state as per bedtools).       

```
# 99460-200 = 99260 (left flank)
# 99460+200 = 99660 (right flank)
# therefore, 99259 - 99660 should capture the 200 bp before and after the variant (in 0-based state)

bedtools getfasta -fi <ref genome> -bed <bed>

# the bed file should be tab delimited, as per
JH816256.1	99259	99660

# this will output to stdout a 401-base record, where the variant is at position 201 (in 1-base state)
>JH816256.1:99259-99660
TACGT....
```
