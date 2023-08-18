# Main genetic analysis script for the pilot study
# Sutherland Bioinformatics
# Initialized 2022-09-18

#### 00. Set up ####
# Source simple_pop_stats and choose Pacific oyster

#### 01. Load Data ####
# Note: if using output of comp_tech_reps, load the RData and your data is obj_nr_best
input.FN <- head(sort(list.files(path = "02_input_data/", pattern = "obj_nr_best", full.names = T), decreasing = T), n = 1)
load(input.FN)

obj_nr_best 
pop(obj_nr_best) # note that all pop attributes are 'unkn' currently, they will be added in the script below
obj <- obj_nr_best

# Otherwise, load the genepop and your data is obj
#load_genepop(datatype = "SNP") 

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

# Use reduced indiv name as indname in genind
indNames(obj) <- indiv.df$indiv

# How many samples from each run? 
table(indiv.df$run)

# Clean-up for write-out
indiv.df <- indiv.df[, "indiv"]
indiv.df <- as.data.frame(indiv.df)
colnames(indiv.df) <- "indiv"
head(indiv.df)

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

## Update population names
# Merge with the population annotation, do not sort
indiv_annot_in_order.df <- merge(x = indiv.df, indiv_annot.df, by = "indiv"
                                 , all.x = T, sort = FALSE # very necessary line
                                 )

head(indiv_annot_in_order.df)

# Observe order remained same as (x) above
head(cbind(indiv_annot_in_order.df, indiv.df), n = 10)
tail(cbind(indiv_annot_in_order.df, indiv.df), n = 10)

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

### TODO: change China to black colour to distinguish

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

abline(h = 50, lty = 3)

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
keep <- obj.df[obj.df$ind.per.missing < 0.5, "ind"]

length(keep)
nInd(obj)

obj.filt <- obj[(keep)]
obj.filt
table(pop(obj.filt))


##### 03.2 Loci - missing data #####
# Filter loci based on missing data
### TODO: this should be generalized into sps ###
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


# Plot
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

# Filter
keep <- rownames(obj.df[obj.df$marker.per.missing < 0.5, ])

# How many loci were removed? 
nLoc(obj.filt)
nLoc(obj.filt) - length(keep)

# Drop loci from genind
obj.all.filt <- obj.filt[, loc=keep]

# Rename back to obj
obj <- obj.all.filt


##### 03.3 Drop monomorphic loci #####
drop_loci(drop_monomorphic = TRUE)
obj <- obj_filt


##### 03.4 Post-QC info collection #####
obj

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

## Density plot stats
# Do high FST markers have lower heterozygosity?
test.df <- per_loc_stats.df
head(test.df)
test.df$elev.fst <- test.df$Fst > 0.05

# Bimodal HOBS
p <- ggplot(test.df, aes(x=Hobs)) +
         geom_density()
p

# Relation between FST and HOBS?
p <- ggplot(test.df, aes(x=Hobs, colour=elev.fst)) +
         geom_density()
p # does not appear to be different
        
            
# Plot Fst by Hobs
pdf(file = "per_locus_Fst_v_Hobs.pdf", width = 8, height = 5)
plot(per_loc_stats.df$Fst, per_loc_stats.df$Hobs
     , las = 1
     , xlab = "Per locus FST"
     , ylab = "Per locus HOBS"
     , pch = 16
     , cex = 0.85
     )

dev.off()



##### 03.5 per marker stats and filters #####
## Per locus statistics
per_locus_stats(data = obj)
head(per_loc_stats.df)

pdf(file = "per_locus_Hobs.pdf", width = 8, height = 5) 
plot(x = per_loc_stats.df$Hobs
     , xlab = "Marker"
     , ylab = "Observed Heterozygosity (Hobs)"
     , las = 1
     , pch = 16
     , cex = 0.85
     )

abline(h = 0.5, lty = 3)
dev.off()

table(per_loc_stats.df$Hobs > 0.6) # only 1

## Optional for dropping Hobs > 0.5 markers ## 
# # Which markers are greater than 0.5 heterozygosity? 
# keep <- setdiff(x = locNames(obj), y = hobs.outliers)
# # Drop Hobs > 0.5 loci from genind
# obj <- obj[, loc=keep]
# obj
## /end/ Optional for dropping ##

## Hardy-Weinberg
hwe_eval(data = obj, alpha = 0.01)


##### 03.5 Post-all filters #####
# Save out colours to be used downstream
colours
colnames(x = colours) <- c("collection", "colour")
write.csv(x = colours, file = "00_archive/formatted_cols.csv", quote = F, row.names = F)


##### 03.6 Relatedness ####
## Identify putative close relatives
# the relatedness script will use the first two values of each individual as the grouping factor, 
# so we need to rename individuals as per their pop IDs
obj.relatedness <- obj
indNames(obj.relatedness) <- paste0(pop(obj), "__", indNames(obj))
indNames(obj.relatedness)

# Calculate inter-individual relatedness
relatedness_calc(data = obj.relatedness, datatype = "SNP") # will output as "03_results/kinship_analysis_<date>.Rdata"
date <- format(Sys.time(), "%Y-%m-%d")
datatype <- "SNP"

# Plot
relatedness_plot(file = paste0("03_results/kinship_analysis_", date, ".Rdata"), same_pops = TRUE, plot_by = "codes", pdf_width = 7, pdf_height = 5)

save.image("03_results/output_relatedness.Rdata")
#load(file = "03_results/output_relatedness.Rdata")

## Inspect relatedness results
gc()

# Load
input.FN <- paste0("03_results/pairwise_relatedness_output_all_", date, ".txt")

# Read in data
rel.df <- read.table(file = input.FN, header = T, sep = "\t")
head(rel.df)

# These are the contrasts in this data
unique(rel.df$group)

# If we ignore the groupings, and just look at the pairs of all individuals, what is the level that outlier is designated?
boxplot(rel.df$ritland, las = 1)
median(rel.df$ritland) # -0.0167
min(boxplot.stats(rel.df$ritland)$out[boxplot.stats(rel.df$ritland)$out > median(rel.df$ritland)]) # 0.1489

## Approach to investigate outliers
# # Set variables of interest
# #popn <- "JPN"
# #popn <- "CHN"
# #popn <- "FRA"
# #popn <- "DPB"
# #popn <- "GUR"
# #popn <- "PEN"
# #popn <- "VIU"
# 
# compare.group <- paste0(substr(x = popn, 1,2), substr(x = popn, 1,2))
# 
# # How many unique inds are there in the selected pop? 
# all_inds.vec  <- c(rel.df$ind1.id, rel.df$ind2.id)
# uniq_inds.vec <- unique(all_inds.vec[grep(pattern = popn, x = all_inds.vec)])
# length(uniq_inds.vec)
# 
# # Select related stat
# statistic <- "ritland"
# 
# # Inspect same-on-same distribution of values
# print(paste0("Identifying outliers using the ", statistic, " statistic"))
# 
# # How many?
# length(rel.df[rel.df$group==compare.group, statistic])
# 
# obs_rel.vec <- rel.df[rel.df$group==compare.group, statistic]
# 
# # Calculate 
# median(obs_rel.vec)
# boxplot.stats(obs_rel.vec)$out
# 
# # Upper outliers
# boxplot.stats(obs_rel.vec)$out[boxplot.stats(obs_rel.vec)$out > median(obs_rel.vec)]
# num_outliers <- length(boxplot.stats(obs_rel.vec)$out[boxplot.stats(obs_rel.vec)$out > median(obs_rel.vec)])
# print(paste0("This relatedness statistic identifies ", num_outliers, " outlier pairs"))
# cutoff <- min(boxplot.stats(obs_rel.vec)$out[boxplot.stats(obs_rel.vec)$out > median(obs_rel.vec)])
# cutoff
# 
# pdf(file = paste0("03_results/related_dist_", compare.group, "_", statistic, ".pdf")
#     , width = 9, height = 5)
# par(mfrow=c(1,2))
# boxplot(obs_rel.vec, las = 1, main = compare.group
#         , ylab = statistic)
# abline(h = cutoff, lty = 3)
# 
# plot(obs_rel.vec, las = 1, main = compare.group
#      , ylab = statistic)
# abline(h = cutoff, lty = 3)
# dev.off()

# Then use the cutoff value to inspect the excel document to ID the pairs, removing one of every two pairs until no outlier pairs remain.

id_close_kin(cutoff = 0.23, statistic = "ritland")
# Will operate on the latest pairwise relatedness oject in the results folder

str(drop.list) # shows how many inds are selected to be dropped from each pop

# Extract all IDs from drop list
drop.inds <- NULL
for(i in 1:length(drop.list)){
  
  drop.inds <- c(drop.inds, drop.list[[i]])
  
}


drop_inds.df <- as.data.frame(drop.inds)
head(drop_inds.df)
drop_inds.df <- separate(data = drop_inds.df, col = "drop.inds", into = c("pop", "ind"), sep = "__", remove = T)
drop.inds <- drop_inds.df$ind
drop.inds

keep.inds <- setdiff(indNames(obj), drop.inds)

obj_purged_relatives <- obj[(keep.inds)]
table(pop(obj_purged_relatives))

table(pop(obj))


#### 04. Analysis ####
## Multivariate
# PCA from genind
pca_from_genind(data = obj_purged_relatives, PCs_ret = 4, colour_file = "00_archive/formatted_cols.csv")

# # DAPC from genind
# dapc_from_genind(data = obj, plot_allele_loadings = TRUE, colour_file = "00_archive/formatted_cols.csv")

## Dendrogram
make_tree(bootstrap = TRUE, boot_obj = obj_purged_relatives, nboots = 10000, dist_metric = "edwards.dist", separated = FALSE)

## Drop low sample size pops
table(pop(obj_purged_relatives))
drop_pops(df = obj_purged_relatives, drop_by_pop_size = TRUE, min_indiv = 15)
obj_pop_filt
table(obj_pop_filt)

## Genetic differentiation
calculate_FST(format = "genind", dat = obj_pop_filt, separated = FALSE, bootstrap = TRUE)


### Polymorphic loci per pop
my_pops.list <- seppop(x = obj, drop=TRUE)

pop_of_interest <- NULL; poly_loci.list <- list(); poly_loci_maf.list <- list()
for(i in 1:length(my_pops.list)){
  
  pop_of_interest <- names(my_pops.list)[i]
  
  print("Dropping non-polymorphic loci")
  
  drop_loci(df = my_pops.list[[i]], drop_monomorphic = TRUE)
  
  poly_loci.list[[pop_of_interest]] <- nLoc(obj_filt)
  
  print("Dropping loci below MAF 0.01")
  maf_filt(data = my_pops.list[[i]], maf = 0.1)
  
  poly_loci_maf.list[[pop_of_interest]] <- nLoc(obj_maf_filt)
  
  
}

poly_loci.list
poly_loci_maf.list


## Private alleles
regional_obj <- obj

# Combine related pops to query private alleles at regional level
unique(pop(regional_obj))
pop(regional_obj) <- gsub(pattern = "VIU_F0|VIU_F1|VIU_F2", replacement = "VIU", x = pop(regional_obj)) # combine VIU
pop(regional_obj) <- gsub(pattern = "PEN|FRA|JPN", replacement = "JPN", x = pop(regional_obj))              # combine JPN lineage
unique(pop(regional_obj))

pa <- private_alleles(gid = regional_obj)
write.csv(x = pa, file = "03_results/private_alleles.csv", quote = F)

save.image(file = "03_results/post-filters_genind_and_enviro.RData")

# Confirm pops are as expected in obj
as.data.frame(cbind(indNames(obj), as.character(pop(obj))))


####### Convert genepop to Rubias format #####
# Need to create a tab-delim stock code file in format of e.g., 
## row 1: collection	repunit
## row 2: boundary_bay	lower_mainland

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
as.data.frame(cbind(indNames(obj), as.character(pop(obj)))) # Note: BR27 should be VIU_F0 [confirmed]
obj # the current analysis object

genepop_to_rubias_SNP(data = obj, sample_type = "reference", custom_format = TRUE, micro_stock_code.FN = micro_stock_code.FN)

print("Your output is available as '03_results/rubias_output_SNP.txt")

file.copy(from = "03_results/rubias_output_SNP.txt", to = "../amplitools/03_prepped_data/cgig_all_rubias.txt", overwrite = T)

save.image(file = "03_results/completed_popgen_analysis.RData")
#load(file = "03_results/completed_popgen_analysis.RData")

# Using this output, move to "amplitools/01_scripts/ckmr_from_rubias.R"


# Additional analysis of the per loc stats related to the reason for selecting markers
perloc.FN <- "./03_results/per_locus_stats_2023-08-16.txt"
marker_reason.FN <- "~/Documents/00_sbio/GBMF_UBC_Pacific_oyster/amplicon_panel/ms_oyster_panel_used_for_design_2022-02-07/04_extract_loci/selected_mnames.csv"

perloc.df <- read.table(file = perloc.FN, header = T, sep = "\t")
head(perloc.df)

reason.df <- read.table(file = marker_reason.FN, header = F, sep = ",")
colnames(reason.df) <- c("mname", "reason")
reason.df <- as.data.frame(reason.df)
head(reason.df)

reason.df <- separate(data = reason.df, col = "mname"
                      , into = c("mname", "pos", "nucl"), sep = "_"
                      , remove = TRUE
                      )
head(reason.df)
reason.df <- reason.df[,c("mname", "reason")]
head(reason.df)

# Combine
all.df <- merge(x = perloc.df, y = reason.df, by = "mname", all.x = T)
head(all.df)

table(all.df[all.df$Hobs > 0.25, "reason"])
length(all.df[all.df$Hobs > 0.25, "reason"])

table(all.df[all.df$Hobs <= 0.25, "reason"])
length(all.df[all.df$Hobs <= 0.25, "reason"])
