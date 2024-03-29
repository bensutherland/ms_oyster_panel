#!/bin/bash
# Uses various .csv files to create a bed file

# Global variables
INPUT_FOLDER="02_input_data"
INPUT_VCF="populations.snps.vcf"
OUTPUT_FOLDER="04_extract_loci"
OUTPUT="vcf_selection.csv"
MNAMES_FN="selected_mnames.csv"

# Remove any existing output
echo "Removing $OUTPUT to make new file" 
rm $OUTPUT_FOLDER/$OUTPUT

# Reporting
echo "*** Duplicate markers reported: ***"
awk -F, '{ print $1 }' $OUTPUT_FOLDER/$MNAMES_FN | 
    sort | uniq -c | sort -n | 
    awk '$1 >1 { print  }' - 

echo "*** Expect the output to have this many markers: ***"
awk -F, '{ print $1 }' $OUTPUT_FOLDER/$MNAMES_FN | 
    sort | uniq | 
    wc -l  

# Prepare bed file
cat $OUTPUT_FOLDER/$MNAMES_FN | 
    
    # Keep only the first column for the mname and remove duplicates
    awk -F, '{ print $1 }' - | 
    sort |
    uniq |

    # Remove positional and variant information from the mname. 
    # Add colon to match VCF. Add mname position, plus 1 due to 0-based count.   
    awk -F_ '{ print $1 ":" $2+1 ":"}' - | 

    # For each mname, extract the relevant section of the vcf
    while read i
    do
        # Reporting
        echo $i

        # Extract the necessary information from the VCF, matching i after the tab
        grep -P "\t"$i $INPUT_FOLDER/$INPUT_VCF |
        awk '{ print $1 "," $2 "," $3 "," $4 "," $5 }' - >> $OUTPUT_FOLDER/$OUTPUT
 
    done    

echo "Completed, the relevant sections of the VCF are in $OUTPUT_FOLDER/$OUTPUT"

