## Additional info on VCF formats 
Some relevant information about formats.      
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

