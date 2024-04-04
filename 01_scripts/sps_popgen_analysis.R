# Main genetic analysis script for the pilot study
# Sutherland Bioinformatics
# Initialized 2022-09-18

#### 00. Set up ####
# Source simple_pop_stats and choose Pacific oyster

# Set variables
max_percent_missing <- 0.3
correct_ids <- TRUE
drop_multimappers <- TRUE
counts_per_locus.FN <- "~/Documents/00_sbio/GBMF_UBC_Pacific_oyster/amplicon_panel/pilot_study/identify_multimappers/counts_per_locus_bowtie2.txt"

#### 01. Load Data ####
# Note: if using output of comp_tech_reps, load the RData and your data is obj_nr_best
input.FN <- head(sort(list.files(path = "02_input_data/", pattern = "obj_nr_best", full.names = T), decreasing = T), n = 1)
load(input.FN)

obj_nr_best 
unique(pop(obj_nr_best)) # all pop attrib currently 'unkn', will be added below
obj <- obj_nr_best
rm(obj_nr_best)

# Otherwise, load the genepop and your data is obj
#load_genepop(datatype = "SNP") 

# the data is here
obj


#### 02. Prepare Data ####
##### 02.1 Manually assign population names based on samples present #####
## Create a list of individuals for manual addition of population
indiv.df <- as.data.frame(indNames(obj))
colnames(indiv.df) <- "indiv"
head(indiv.df)

# Separate components of indiv ID (i.e., run, barcode, sample)
indiv.df <- separate(data = indiv.df, col = "indiv", into = c("run", "barcode", "indiv"), sep = "__", remove = T)
head(indiv.df)

# How many samples from each run? 
table(indiv.df$run)

# Use reduced indiv name as indname in genind
indNames(obj) <- indiv.df$indiv

##### Update individual names #####
# Based on an earlier analysis of parentage, several individuals were determined to be 
#  mislabeled, or were from a mixed family. Therefore we will correct these individuals here
#  based on the variable above 'correct_ids=TRUE'
if(correct_ids==TRUE){
  
  # Correct the known erroneous samples as follows
  inds <- indNames(obj)
  inds <- gsub(pattern = "F10-01", replacement = "F9-101", inds)
  inds <- gsub(pattern = "F14-30", replacement = "F17-130", inds)
  inds <- gsub(pattern = "F2-11", replacement = "F3-111", inds)
  inds <- gsub(pattern = "F15-43", replacement = "F14-143", inds)
  inds <- gsub(pattern = "F15-10", replacement = "F14-110", inds)
  inds <- gsub(pattern = "F15-06", replacement = "F14-106", inds)
  
  # Update the indiv names in the genind
  indNames(obj) <- inds
  
}else{
  
  print("Not changing IDs. Running analysis with known erroneous sample labels.")
  
}

## Create a list of individuals for manual addition of population
indiv.df <- as.data.frame(indNames(obj))
colnames(indiv.df) <- "indiv"
head(indiv.df)

# Clean-up for write-out
indiv.df <- as.data.frame(indiv.df)
head(indiv.df)
dim(indiv.df)

# Add dummy column to fill manually
indiv.df$pop <- NA

# Write out empty file to provide pop names
write.table(x = indiv.df, file = "02_input_data/my_data_ind-to-pop.txt"
            , sep = "\t", col.names = T, row.names = F
            , quote = F
            )

# In folder above, *manually annotate* the output file above
# , save with "_annot.txt" appended, populate with pop names (no spaces)

# Load annotated df
indiv_annot.df <- read.table(file = "02_input_data/my_data_ind-to-pop_annot.txt"
                           , header = T, sep = "\t"
                           #, quote = F
                           )
dim(indiv_annot.df)
head(indiv_annot.df)

## Update population names
# Merge with the population annotation, do not sort
indiv_annot_in_order.df <- merge(x = indiv.df, indiv_annot.df, by = "indiv"
                                 , all.x = T, sort = FALSE # very necessary line
                                 )

head(indiv_annot_in_order.df)
dim(indiv_annot.df)

# Observe order remained same as (x) above
head(cbind(indiv_annot_in_order.df, indiv.df), n = 20)
tail(cbind(indiv_annot_in_order.df, indiv.df), n = 20)

### TODO: add data-check in this step###
# # Write a little test to be sure
# test <- cbind(indNames(obj), indiv_annot_in_order.df$indiv)
# table(test[,1] == test[ ,2])
# 
# # test <- cbind(indNames(obj), sort(indiv_annot_in_order.df$indiv))
# # table(test[,1] == test[ ,2])
# ## /END/ ##

# Assign the pop IDs to the genind
pop(obj) <- indiv_annot_in_order.df$pop.y
table((pop(obj)))

##### 02.2 Add in population colours #####
## Population colours
pops_in_genepop <- unique(pop(obj))
pops_in_genepop.df <- as.data.frame(pops_in_genepop)

## Download colours file from previous git repo
# url = "https://raw.githubusercontent.com/bensutherland/ms_oyster_popgen/master/00_archive/my_cols.csv" # only need to run once
destfile <- "00_archive/my_cols.csv"
# download.file(url, destfile)   # only need to run once
my_colours <- read.csv(destfile)
new_pop_colours <- matrix(c("VIU_F0", "VIU_F1", "VIU_F2", "darkred", "red", "magenta"), nrow = 3, ncol = 2)
colnames(new_pop_colours) <- c("my.pops", "my.cols")
my_colours <- rbind(my_colours, new_pop_colours)
my_colours

# Connect colours to empirical populations, which will exclude any that are not in the dataset
colours <- merge(x = pops_in_genepop.df, y =  my_colours, by.x = "pops_in_genepop", by.y = "my.pops"
                 #, sort = F
                 , all.x = T
                 )
colours

# Clean space
rm(my_colours)
rm(new_pop_colours)


#### 03. Characterize missing data (indiv and loci) and filter ####
# Set variables to use for both plots
plot_width  <- 8
plot_height <- 5
plot_cex    <- 0.85
plot_pch    <- 16

##### 03.1 Individuals - missing data #####
percent_missing_by_ind(df = obj)
head(missing_data.df)

# Add pop IDs to the missing data df
missing_data.df <- merge(x = missing_data.df, y = indiv_annot.df, by.x = "ind", by.y = "indiv", all.x = T)
head(missing_data.df)
head(indiv_annot.df)

# Add colours to the missing data df
colours
plot_cols.df <- merge(x = missing_data.df, y = colours, by.x = "pop", by.y = "pops_in_genepop", all.x = T
                      , sort = F
                      )

# Plot missing data by individual
pdf(file = "03_results/geno_rate_by_ind.pdf", width = plot_width, height = plot_height)
plot(100 * (1 - plot_cols.df$ind.per.missing), ylab = "Genotyping rate (%)"
     , col = plot_cols.df$my.cols
     , las = 1
     , xlab = "Individual"
     , ylim = c(0,100)
     , pch=plot_pch
     , cex = plot_cex
     )

abline(h = 70, lty = 3)

legend("bottomleft", legend = unique(plot_cols.df$pop)
       , fill = unique(plot_cols.df$my.cols)
       , cex = 0.7
       , bg = "white"
       )
dev.off()

# Keep missing data info on retained indiv
obj.df <- missing_data.df
head(obj.df)

## Filter individuals by genotyping rate 
keep <- obj.df[obj.df$ind.per.missing < max_percent_missing, "ind"]

length(keep)
nInd(obj)

obj.filt <- obj[(keep)]
obj.filt
table(pop(obj.filt))


##### 03.2 Loci - missing data #####
# Filter loci based on missing data
obj.df <- genind2df(obj.filt)
obj.df[1:5,1:5]
obj.df <- t(obj.df)
obj.df[1:5,1:5]
obj.df <- obj.df[2:nrow(obj.df),] # remove pop row
obj.df[1:5,1:5]
dim(obj.df)
str(obj.df)

obj.df <- as.data.frame(obj.df)
dim(obj.df)
#str(obj.df)
obj.df[1:5,1:5] # See top left of file
obj.df[(dim(obj.df)[1]-5):dim(obj.df)[1], (dim(obj.df)[2]-5):dim(obj.df)[2]] # See bottom right of file

# Add collector col
obj.df$marker.per.missing <- NA

# Loop across df, checking each row for the number of NAs and dividing by the total number of markers in the df
for(i in 1:(nrow(obj.df))){
  
  # Per marker                      sum all NAs for the marker, divide by total number markers
  obj.df$marker.per.missing[i] <-  (sum(is.na(obj.df[i,]))-1) / (ncol(obj.df)-1) 
  
}


# Plot
pdf(file = "03_results/geno_rate_by_marker.pdf", width = plot_width, height = plot_height)
plot(100 * (1- obj.df$marker.per.missing), xlab = "Marker", ylab = "Genotyping rate (%)", las = 1
     , ylim = c(0,100)
     , pch = plot_pch
     , cex = plot_cex
     )
abline(h = 70
       #, col = "grey60"
       , lty = 3)
dev.off()

# Filter
keep <- rownames(obj.df[obj.df$marker.per.missing < max_percent_missing, ])

# How many loci were removed? 
nLoc(obj.filt)
print(paste0("Dropping ", (nLoc(obj.filt) - length(keep)), " markers"))

# Drop loci from genind
obj.all.filt <- obj.filt[, loc=keep]

# Rename back to obj
obj <- obj.all.filt
rm(obj.all.filt)
obj


##### 03.3 Drop monomorphic loci #####
drop_loci(drop_monomorphic = TRUE)
obj <- obj_filt
rm(obj_filt)
rm(obj.filt)
gc()
obj


##### 03.4 Drop multi-mappers #####
if(drop_multimappers==TRUE){
  
  # Read in the number alignments per locus file
  counts_per_locus.FN
  counts_per_locus.df <- read.table(file = counts_per_locus.FN, header = F)
  head(counts_per_locus.df)
  colnames(counts_per_locus.df) <- c("count", "mname")
  
  # How many have more than one alignment? 
  loci_to_remove.vec <- counts_per_locus.df[counts_per_locus.df$count > 1, "mname"]
  length(loci_to_remove.vec)
  write.table(x = loci_to_remove.vec, file = "03_results/multimapper_loci_to_remove.txt", quote = F, sep = "\t" , row.names = F, col.names = F)
  
  drop_loci(df = obj, drop_monomorphic = F, drop_file = "03_results/multimapper_loci_to_remove.txt")
  obj <- obj_filt
  
  
}else{
  
  print("Note: not removing multi-mappers")
  
}
obj


##### 03.5 per marker stats and filters #####
## Per locus statistics
per_locus_stats(data = obj)
head(per_loc_stats.df)
nrow(per_loc_stats.df)

## Density plot stats
# Do high FST markers have lower heterozygosity?
head(per_loc_stats.df)
#table(per_loc_stats.df$Fst > 0.05)

# Plot Fst by Hobs
pdf(file = "03_results/per_locus_Fst_v_Hobs.pdf", width = 8, height = 5)
plot(per_loc_stats.df$Fst, per_loc_stats.df$Hobs
     , las = 1
     , xlab = "Per locus FST"
     , ylab = "Per locus HOBS"
     , pch = 16
     , cex = 0.85
)

dev.off()

# How many loci have more than 0.5 HOBS? 
table(per_loc_stats.df$Hobs > 0.5)

## Hardy-Weinberg, only on wild samples though ##
# Remove cultured CHN samples, previously unknown but now noted as QDC samples
QDC.inds <- c("1601", "1602", "1603", "1604", "1605", "1606", "1607", "1608")
keep.inds <- setdiff(x = indNames(obj), y = QDC.inds)
obj.no.QDC <- obj[keep.inds]

hwe_eval(data = obj.no.QDC, alpha = 0.01)

# Identify hwe outliers in the 'wild' pops
col.oi <- grep(pattern = "Pr", x = colnames(per_locus_hwe_JPN.df))
hwe_outlier_mname_JPN.vec <- per_locus_hwe_JPN.df[per_locus_hwe_JPN.df[, col.oi] < 0.01, "mname"]
hwe_outlier_mname_PEN.vec <- per_locus_hwe_PEN.df[per_locus_hwe_PEN.df[, col.oi] < 0.01, "mname"]
hwe_outlier_mname_FRA.vec <- per_locus_hwe_FRA.df[per_locus_hwe_FRA.df[, col.oi] < 0.01, "mname"]
hwe_outlier_mname_CHN.vec <- per_locus_hwe_CHN.df[per_locus_hwe_CHN.df[, col.oi] < 0.01, "mname"]

length(hwe_outlier_mname_JPN.vec)
length(hwe_outlier_mname_PEN.vec)
length(hwe_outlier_mname_FRA.vec)
length(hwe_outlier_mname_CHN.vec)

# Collect the names of the HW deviators
hwe_outlier.df <- as.data.frame(sort(table(c(hwe_outlier_mname_JPN.vec, hwe_outlier_mname_PEN.vec, hwe_outlier_mname_FRA.vec, hwe_outlier_mname_CHN.vec)), decreasing = T))
colnames(hwe_outlier.df) <- c("mname", "freq")
head(hwe_outlier.df, n = 20)
nrow(hwe_outlier.df)
write.table(x = hwe_outlier.df, file = "03_results/hwe_outlier_summary.txt", quote = F, sep = "\t", row.names = F)


##### 03.6 Post-QC info collection #####
obj


## View the ind or loc names
inds <- indNames(obj)
loci <- locNames(obj)

# Save out which individuals have passed the filters
write.table(x = inds, file = "03_results/post-popgen_filters_retain_indiv.txt", sep = "\t", quote = F
            , row.names = F, col.names = F
            )

write.table(x = loci, file = "03_results/post-popgen_filters_retain_loci.txt", sep = "\t", quote = F
            , row.names = F, col.names = F
)



##### 03.7 Post-all filters #####
# Save out colours to be used downstream
colours
colnames(x = colours) <- c("collection", "colour")
write.csv(x = colours, file = "00_archive/formatted_cols.csv", quote = F, row.names = F)

# Preserve all data before moving into parentage analysis
obj.bck <- obj


# Save output
save.image("03_results/post-filters_popgen_dataset.RData")



# Or continue below


#### 04. Parentage: convert genepop to Rubias format ####
# Limit the data to only include the VIU pops
separated_pops <- seppop(obj)
obj_parentage <- repool(separated_pops$VIU_F0, separated_pops$VIU_F1, separated_pops$VIU_F2)
table(pop(obj_parentage))
obj <- obj_parentage

## Aside, instruction: need to create a tab-delim stock code file in following format
# row 1: # collection	repunit
# row 2: # boundary_bay	lower_mainland

# Here we will just create a df based on existing populations where collection = repunit
stock_code.df <- as.data.frame(unique(pop(obj)))
colnames(stock_code.df) <- "collection"
stock_code.df$repunit <- stock_code.df$collection
stock_code.df

# Write it out
write_delim(x = stock_code.df, file = "00_archive/stock_code.txt", delim = "\t", col_names = T)
micro_stock_code.FN <- "00_archive/stock_code.txt"
# this is for annotate_rubias(), for an unknown reason it requires the name micro_stock_code.FN

## Convert genepop to rubias
datatype <- "SNP" # required for genepop_to_rubias_SNP
#as.data.frame(cbind(indNames(obj), as.character(pop(obj))))
#obj # the current analysis object

## genepop to rubias with all loci, including multi-mappers
# genepop_to_rubias_SNP(data = obj, sample_type = "reference", custom_format = TRUE, micro_stock_code.FN = micro_stock_code.FN
#                       , pop_map.FN = "02_input_data/my_data_ind-to-pop_annot.txt")
# 
# print("Your output is available as '03_results/rubias_output_SNP.txt")
# 
# file.copy(from = "03_results/rubias_output_SNP.txt", to = "../amplitools/03_prepped_data/cgig_all_rubias.txt", overwrite = T)

# Remove monomorphic sites
obj
drop_loci(df = obj, drop_monomorphic = T)
obj <- obj_filt


#### Retain information about loci kept for parentage ####
## View the ind or loc names
inds <- indNames(obj)
length(inds)
loci <- locNames(obj)
length(loci)

# Save out which individuals have passed the filters
write.table(x = inds, file = "03_results/post-filters_parentage_retain_indiv.txt", sep = "\t", quote = F
            , row.names = F, col.names = F
)

write.table(x = loci, file = "03_results/post-filters_parentage_retain_loci.txt", sep = "\t", quote = F
            , row.names = F, col.names = F
)


## All filtered loci: write to rubias for parentage assignment
pop_map.FN <- "02_input_data/my_data_ind-to-pop_annot.txt"
genepop_to_rubias_SNP(data = obj, sample_type = "reference"
                      , custom_format = TRUE, micro_stock_code.FN = micro_stock_code.FN
                      , pop_map.FN = pop_map.FN
)
print("Your output is available as '03_results/rubias_output_SNP.txt")

# Obtain some variables to create filename
date <- format(Sys.time(), "%Y-%m-%d")
indiv.n <- nInd(obj)
loci.n  <- nLoc(obj)
rubias_custom.FN <- paste0("03_results/rubias_", indiv.n, "_ind_", loci.n, "_loc_", date, ".txt")

# Copy to rename rubias output file
file.copy(from = "03_results/rubias_output_SNP.txt", to = rubias_custom.FN, overwrite = T)
print(paste0("The output is saved as ", rubias_custom.FN))

# Save image
save.image("03_results/post-filters_prepared_for_parentage_rubias_built.RData")

print("Here you need to copy the above rubias file to amplitools results folder.")

# Using this output, move to "01_scripts/sps_popgen_analysis_part_2_parentage.R"
# OR
# move to "01_scripts/sps_popgen_analysis_part_2_relatedness.R" for additional popgen analyses
# to analyze the relatives and post-purged sibs analysis in the script
