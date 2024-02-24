## Analysis script for OCP23 (parentage analysis)
# Ben Sutherland and Liam Surry, VIU
# initialized 2023-08-24

# If you are running from raw data, do step 1. If you've already prepared a genepop, go to step 2. 

#### 01. Prepare genepop in amplitools ####
# Open and source amplitools/01_scripts/00_initiator.R

# Convert proton to genepop
proton_to_genepop(neg_control = "Blank")

# Prepare the genepop
#./01_scripts/format_genepop.sh 02_input_data/prepped_matrices/<filename>.txt

# Copy the results to simple_pop_stats
# cp ./02_input_data/prepped_genepops/R_2023_07_26_12_44_23_user_GSS5PR-0268-78-Ampseq_Oyster_20230725_gen_data.gen ../simple_pop_stats/02_input_data/


#### 02. Filter in simple_pop_stats ####
# Open and source simple_pop_stats/01_scripts/simple_pop_stats_start.R, then choose '9-Pacific oyster'
## note: change variable 'on_network' to FALSE

#### 01. Load data ####
load_genepop(datatype = "SNP")
## note: input file is "02_input_data/R_2023_07_26_12_44_23_user_GSS5PR-0268-78-Ampseq_Oyster_20230725_gen_data.gen"
head(indNames(obj)) # indiv names are in standard amplitools format

#### 02. Prepare data ####
# Simplify amplitools names
simplify_names(df = obj, format = "amplitools")
obj <- obj_simplified
indNames(obj)


##### 02.1 Manually assign population names based on samples present #####
generate_popmap(df = obj)

# Manually annotate the pop map file "02_input_data/my_data_ind-to-pop.txt"
# , save with "_annot.txt" appended, populate with pop names (no spaces)

## Load annotated df
indiv_annot.df <- read.table(file = "00_archive/my_data_ind-to-pop_annot.txt"
                             , header = T, sep = "\t"
                             #, quote = F
)
head(indiv_annot.df)

## Update pop attribute
### TODO: this should be a function ##
# Obtain sample names from obj and keep as indiv.df
indiv.df <- NULL
indiv.df <- indNames(obj)
indiv.df <- as.data.frame(indiv.df)
colnames(indiv.df) <- "indiv"

# These are the two files that will be merged
head(indiv.df)
head(indiv_annot.df)

# Merge the ordered sample names with updated population annotation, do not sort
indiv_annot_in_order.df <- merge(x = indiv.df, indiv_annot.df, by = "indiv"
                                 , all.x = T, sort = FALSE # very necessary line
)

head(indiv_annot_in_order.df)
tail(indiv_annot_in_order.df)

# Observe order remained same as (x) above
head(cbind(indiv_annot_in_order.df, indiv.df), n = 10)
tail(cbind(indiv_annot_in_order.df, indiv.df), n = 10)

# Update the pop attribute from the ordered sample metadata
pop(obj) <- indiv_annot_in_order.df[, "pop"]
table((pop(obj)))


##### 02.2 Set population colours #####
## Population colours
colours <- matrix(c("F1", "F0", "OAR", "darkgreen", "purple", "black"), nrow = 3, ncol = 2)
colnames(colours) <- c("my.pops", "my.cols")
colours

# Save out colours to be used downstream
colnames(x = colours) <- c("collection", "colour")
write.csv(x = colours, file = "00_archive/formatted_cols.csv", quote = F, row.names = F)


#### 03. Characterize missing data (indiv and loci) and filter ####
### This should be a function (SPS)
# Set variables to use for both plots
plot_width <- 8
plot_height <- 5
plot_cex <- 0.85
plot_pch <- 16

##### 03.1 Individuals - missing data #####
# Set variables
max_missing <- 0.3

percent_missing_by_ind(df = obj)
head(missing_data.df)

# Add pop attribute to the missing data df
head(missing_data.df)
head(indiv_annot.df)
missing_data.df <- merge(x = missing_data.df, y = indiv_annot.df, by.x = "ind", by.y = "indiv", all.x = T)
head(missing_data.df)
dim(missing_data.df)

# Can observe missing data by sex
boxplot(missing_data.df$ind.num.typed ~ missing_data.df$sex)

# Add colours to the missing data df
colours
plot_cols.df <- merge(x = missing_data.df, y = colours, by.x = "pop", by.y = "collection", all.x = T
                      , sort = F
)
head(plot_cols.df)

# Plot missing data by individual
pdf(file = "03_results/geno_rate_by_ind.pdf", width = plot_width, height = plot_height)
plot(100 * (1 - plot_cols.df$ind.per.missing), ylab = "Genotyping rate (%)"
     , col = plot_cols.df$colour
     , las = 1
     , xlab = "Individual"
     , ylim = c(0,100)
     , pch=plot_pch
     , cex = plot_cex
)

abline(h = 50, lty = 3)

legend("bottomleft", legend = unique(plot_cols.df$pop)
       , fill = unique(plot_cols.df$colour)
       , cex = 0.7
       , bg = "white"
)
dev.off()

# Filter individuals by missing data
obj.df <- missing_data.df
head(obj.df)
dim(obj.df)

## Remove duplicate individuals
# sort by percent missing
obj.df <- obj.df[with(obj.df, order(ind.per.missing)), ]
head(obj.df, n = 20)
obj.df <- obj.df[!duplicated(obj.df$alt.ID),] # removes the duplicate with most missing
dim(obj.df)
table(obj.df$pop)

keep <- obj.df[obj.df$ind.per.missing < max_missing, "ind"]

length(keep)
nInd(obj)

obj.filt <- obj[(keep)]
obj.filt

# Samples remaining after filters
table(pop(obj.filt))

## Remove unrelated project samples
obj.sep <- seppop(x = obj.filt)
obj.filt <- repool(obj.sep$F1, obj.sep$F0)
table(pop(obj.filt))
obj.filt


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
str(obj.df)
obj.df[1:5,1:5] # See top left of file
obj.df[(dim(obj.df)[1]-5):dim(obj.df)[1], (dim(obj.df)[2]-5):dim(obj.df)[2]] # See bottom right of file

# Add collector col
obj.df$marker.per.missing <- NA

for(i in 1:(nrow(obj.df))){
  
  # Per marker                      sum all NAs for the marker, divide by total number markers
  obj.df$marker.per.missing[i] <-  (sum(is.na(obj.df[i,]))-1) / (ncol(obj.df)-1) 
  
}


# Plot missing data by marker
pdf(file = "03_results/geno_rate_by_marker.pdf", width = plot_width, height = plot_height)
plot(100 * (1- obj.df$marker.per.missing), xlab = "Marker", ylab = "Genotyping rate (%)", las = 1
     , ylim = c(0,100)
     , pch = plot_pch
     , cex = plot_cex
)
abline(h = 50
       #, col = "grey60"
       , lty = 3)
dev.off()

# Filter markers by genotyping rate
keep <- rownames(obj.df[obj.df$marker.per.missing < max_missing, ])

# How many loci will be removed? 
nLoc(obj.filt)
nLoc(obj.filt) - length(keep)

# Drop loci from genind
obj.all.filt <- obj.filt[, loc=keep]

# Rename back to obj
obj <- obj.all.filt
obj


##### 03.3 Drop monomorphic loci #####
drop_loci(df = obj, drop_monomorphic = TRUE) # drops monomorphic markers
obj <- obj_filt

# Final remaining dataset
obj


##### 03.5 per marker stats #####
# MAF information
maf_filt(data = obj, maf = 0.0001) # for generating myFreq to plot only
obj <- obj_maf_filt
head(myFreq)

pdf(file = "03_results/maf_freq.pdf", width = 7, height = 5)
hist(myFreq, breaks = 20, main = "", xlab = "MAF", las = 1)
dev.off()

## Per locus statistics
per_locus_stats(data = obj)
head(per_loc_stats.df)

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

## View the ind or loc names
inds <- indNames(obj)
loci <- locNames(obj)

# Save out which individuals have passed the filters
write.table(x = inds, file = "03_results/retained_individuals.txt", sep = "\t", quote = F
            , row.names = F, col.names = F
)

write.table(x = loci, file = "03_results/retained_loci.txt", sep = "\t", quote = F
            , row.names = F, col.names = F
)


# Summary of excess HOBS
table(per_loc_stats.df$Hobs > 0.5) 
# note: due to family nature of the data, not filtering based on HWP or HOBS 


##### 03.6 Drop multi-mappers #####
counts_per_locus.FN <- "~/Documents/00_sbio/GBMF_UBC_Pacific_oyster/amplicon_panel/pilot_study/identify_multimappers/counts_per_locus.txt"
counts_per_locus.df <- read.table(file = counts_per_locus.FN, header = F)
head(counts_per_locus.df)
colnames(counts_per_locus.df) <- c("count", "mname")

# How many have more than one? 
loci_to_remove.vec <- counts_per_locus.df[counts_per_locus.df$count > 1, "mname"]
length(loci_to_remove.vec)
write.table(x = loci_to_remove.vec, file = "03_results/loci_to_remove.txt", quote = F, sep = "\t" , row.names = F, col.names = F)

drop_loci(df = obj, drop_monomorphic = F, drop_file = "03_results/loci_to_remove.txt")

# Save output
save.image(file = "03_results/filtered_genind_before_ckmr.RData")

# Next, go to OCP23_analysis_part_2_2024-02-23.R 
