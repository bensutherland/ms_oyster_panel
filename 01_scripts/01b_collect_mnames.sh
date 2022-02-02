#!/bin/bash
# Uses various .csv files to collect mnames 

# Global variables
INPUT_FOLDER="03_marker_selection"
INPUT_FN="all_markers_per_locus_stats_incl_hobs_MAF.csv"
OUTPUT_FOLDER="04_extract_loci"
OUTPUT_FN="selected_mnames.csv"
NUM_MARKERS=300

# Remove any existing output
echo "Removing $OUTPUT_FN to make new file" 
rm $OUTPUT_FOLDER/$OUTPUT_FN

# Obtain FST markers
echo "Taking the $NUM_MARKERS highest FST markers"
grep -vE '^mname' $INPUT_FOLDER/$INPUT_FN | 
    awk -F, '{ print $1 "," $3 }' - | 
    sort -t, -gk2 -r - |
    awk -F, '{ print $1 ", top_FST" }' - | 
    head -n $NUM_MARKERS > $OUTPUT_FOLDER/$OUTPUT_FN  

# Obtain Hobs markers
echo "Taking the $NUM_MARKERS highest Hobs (under 0.5) markers"
grep -vE '^mname' $INPUT_FOLDER/$INPUT_FN | 
    awk -F, '$6<=0.5 { print $1 "," $6 }' - | 
    sort -t, -gk2 -r - | 
    awk -F, '{ print $1 ", top Hobs" }' - |
    head -n $NUM_MARKERS >> $OUTPUT_FOLDER/$OUTPUT_FN

# Obtain custom markers
echo "Taking custom mnames and adding to the collection"
cat $INPUT_FOLDER/*selected*mnames.csv |
    awk -F, '{ print $1 ", CUSTOM"}' - >> $OUTPUT_FOLDER/$OUTPUT_FN

echo "Finished, selected mnames are in $OUTPUT_FOLDER/$OUTPUT_FN"


