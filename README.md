# ms_oyster_panel
Repository for code to support the manuscript to design a Pacific oyster amplicon panel

Note: this repository is for the designed use of the author only and comes with no guarantees of usefulness for anything else.       

Requirements:      
enormandeau/Scripts     

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

Marker file: 
```
mname,Fit,Fst,Fis,maf.vec,Hobs
100026_37_C,-0.011,0.00777,-0.0192,0.0173,0.0347 
```


### 02. Define the markers to be used for the extraction ###
Use the marker file from above, and select top FST or top Hobs then extract the relevant information from the VCF.      
Sort the marker file, then select the top 300 FST markers:     
`grep -vE '^mname' 03_marker_selection/all_markers_per_locus_stats_incl_hobs_MAF.csv | awk -F, '{ print $1 "," $3 }' - | sort -t, -gk2 -r - | head -n 300 > 04_extract_loci/top_fst_mname_fst.csv`        

Output will be rows of markers, no header, e.g., `585541_17_C,0.082`      


Sort the marker file, keep only columns where Hobs is less than or equal 0.5, then select the top 300 Hobs markers:     
`awk -F, '$6<=0.5 { print $1 "," $6 }' 03_marker_selection/all_markers_per_locus_stats_incl_hobs_MAF.csv | sort -t, -gk2 -r - | head -n 300 > 04_extract_loci/top_hobs_mname_hobs.csv`      

Output will be rows of markers, no header, e.g., `92947_23_A,0.5`     


### 03. Use the defined markers to extract the relevant lines from the VCF
The above output will be used to extract the relevant lines and fields from the VCF. To conduct this, use the following script:        
`01_scripts/02_info_from_vcf.sh`          
...this will output `04_extract_loci/vcf_selection.csv`        

e.g.,     
```
JH816256.1,99460,100388:26:-,G,A

# where the fields are: 
# (1) chromosome
# (2) position of the SNP in the reference genome
# (3) information about the marker
# (4) reference allele
# (5) variant allele

```


To extract from the reference genome, we will use bedtools combined with a bedfile.     
The bedfile will be in the shape of:     
```
chr1 5 10
```
...where the file is space separated, and the positions refer to positions in the contig/scaffold starting counting from position 0.      


Are we sure about what we are looking at?      

To get a better view of the fasta file, use the enormandeau/Scripts to unwrap the fasta file so that each record (post fasta header) takes up a single line.     
`fasta_unwrap.py 02_input_data/GCA_000297895.1_oyster_v9_genomic.fna 02_input_data/GCA_000297895.1_oyster_v9_genomic_unwrapped.fna`     

Pull out an example
`JH816256.1,99460,100388:26:-,G,A` from `04_extract_loci/vcf_selection.csv`       

Let's take the fasta consensus sequence for this from the stacks fasta file as well, for extra information.       
```
>CLocus_100388 [JH816256.1, 99396, -]
TGCATCTGTCCTCTTTTCTGTGCTTCATTGATGGTCATTTTGGATACCTGTTTGACTTAACTGATTTAGATAAGATGATCATGTGTTGTG
```

From the genome, it looks like:     
`grep -A1 JH816256.1 02_input_data/GCA_000297895.1_oyster_v9_genomic_unwrapped.fna > 02_input_data/example_record_JH816256.1.txt`     

Looking at the reference genome record (open in text editor) you will see that the position 99396 (stacks fasta) indicates the starting character of the reverse complement of the stacks fasta record.     
Further, looking at position 99460 in the ref genome, you will see that this is the 'G' nucleotide, which is the reference allele (i.e., TGAAGC). Note: this would be in a 1-based orientation.        

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

Note that the above will not insert square brackets, so this will need to be done later.    

# 04 Extract FASTA from reference genome
Use the following to prepare a bed file from the relevant lines of the bcf
`01_scripts/03_prepare_bed_file.sh`     


