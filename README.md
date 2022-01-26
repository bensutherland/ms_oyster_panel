# ms_oyster_panel
Code to support the design of a Pacific oyster amplicon panel.     

Note: this repository is for the designed use of the author only and comes with no guarantees of usefulness for anything else.       

Requirements:      
enormandeau/Scripts     

Data inputs:     
- Single SNP data in plink format from Sutherland et al. 2020 (Evol. Appl.)        
- Single SNP data populations output in VCF format       
- Reference genome for Pacific oyster used in identifying markers (Zhang et al. 2012; Sutherland et al. 2020)         
Download this version, GenBank format: https://www.ncbi.nlm.nih.gov/assembly/GCF_000297895.1/        

Note: in Sutherland et al. 2020, the genetic data was been filtered to retain markers present in 70% of individuals in all of 15 populations, with a --min-maf >= 0.01, and with HWE outliers removed.        

All shell scripts are run from the main repo.     


### 01. Filter and characterize variants ###
Use the script `01_scripts/01_identifying_markers.R` run interactively in R. This script will do the following:      
1. Load packages and set working directory
2. Import data and select only the BC naturalized samples
3. Perform minor allele frequency (MAF) filter to remove low MAF variants wrt. BC pops
4. Calculate per locus stats observed heterozygosity (Hobs) and FST (averaged per locus)
5. Generate comparative plots
6. Output a csv of marker per locus stats for all MAF-filtered variants as `03_marker_selection/all_markers_per_locus_stats_incl_hobs_MAF.csv`     

Marker file: 
```
mname,Fit,Fst,Fis,maf.vec,Hobs
100026_37_C,-0.011,0.00777,-0.0192,0.0173,0.0347 
```


### 02. Select the top markers to create whitelists ###
Use the marker file from above, and select top FST or top Hobs then extract the relevant information from the VCF.      
Top FST:     
`grep -vE '^mname' 03_marker_selection/all_markers_per_locus_stats_incl_hobs_MAF.csv | awk -F, '{ print $1 "," $3 }' - | sort -t, -gk2 -r - | head -n 300 > 04_extract_loci/top_fst_mname_fst.csv`        

Output will be rows of markers, no header, e.g., `585541_17_C,0.082`      

Top Hobs (but below 0.5):        
`awk -F, '$6<=0.5 { print $1 "," $6 }' 03_marker_selection/all_markers_per_locus_stats_incl_hobs_MAF.csv | sort -t, -gk2 -r - | head -n 300 > 04_extract_loci/top_hobs_mname_hobs.csv`      

Output will be rows of markers, no header, e.g., `92947_23_A,0.5`     


### 03. Collect the whitelist marker info from VCF ###
Run the following script:       
`01_scripts/02_info_from_vcf.sh`          
...this will output `04_extract_loci/vcf_selection.csv`, which contains (for example):             


```
JH816256.1,99460,100388:26:-,G,A

# where the fields are: 
# (1) chromosome
# (2) position of the SNP in the reference genome
# (3) information about the marker
# (4) reference allele
# (5) variant allele

```

### 04. Overview ####
This section is optional, and shows the logic underlying following steps.      
To extract from the reference genome, we will use bedtools combined with a bedfile.     
The bedfile will be in the shape of:     
```
chr1	5	10
```
...where the file is tab-delimited, and the positions refer to selecting positions in the contig/scaffold starting counting from position 0.      

To get a better view of the fasta file, use the enormandeau/Scripts to unwrap the fasta file so that each record's nucleotide section takes up a single line.     
`fasta_unwrap.py 02_input_data/GCA_000297895.1_oyster_v9_genomic.fna 02_input_data/GCA_000297895.1_oyster_v9_genomic_unwrapped.fna`     

Take an example:       
`JH816256.1,99460,100388:26:-,G,A` from `04_extract_loci/vcf_selection.csv`       

Also review the fasta consensus sequence for this from the stacks fasta file:      
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


# 04 Extract FASTA from reference genome
Use the following to prepare a bed file from the relevant lines of the bcf
`01_scripts/03_prepare_bed_file.sh`     


