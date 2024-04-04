# Relatedness analysis script and post-purged sibs analyses for the pilot study
# Sutherland Bioinformatics
# Initialized 2022-09-18

#### 00. Set up ####
# Source simple_pop_stats and choose Pacific oyster

# User set variables
marker_reason.FN <- "~/Documents/00_sbio/GBMF_UBC_Pacific_oyster/amplicon_panel/pilot_study/z-draft_analyses/ms_oyster_panel_used_for_design_2022-02-07/04_extract_loci/selected_mnames.csv"

#load(file = "03_results/output_post-filters.Rdata")
load(file = "03_results/post-filters_prepared_for_parentage_rubias_built.RData")

# Working with the following: 
obj.bck # 312 inds, 409 loci
table(pop(obj.bck)) # this has all inds

obj <- obj.bck

##### Private alleles ####
regional_obj <- obj

# Combine similar collections prior to identifying private alleles
unique(pop(regional_obj))
pop(regional_obj) <- gsub(pattern = "VIU_F0|VIU_F1|VIU_F2", replacement = "VIU", x = pop(regional_obj)) # combine VIU
pop(regional_obj) <- gsub(pattern = "PEN|FRA|JPN", replacement = "JPN", x = pop(regional_obj))              # combine JPN lineage
unique(pop(regional_obj))
table(pop(regional_obj))

pa <- private_alleles(gid = regional_obj)
write.csv(x = pa, file = "03_results/private_alleles.csv", quote = F)


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
median(rel.df$ritland) # -0.0165
min(boxplot.stats(rel.df$ritland)$out[boxplot.stats(rel.df$ritland)$out > median(rel.df$ritland)]) # 0.1533

# Check full-sibs of VIU pops
col1 <- grep(pattern = "VIU_F2", x = rel.df$ind1.id)
col2 <- grep(pattern = "VIU_F2", x = rel.df$ind2.id)
keep_rows <- c(col1, col2)
rel_VIU_F2.df <- rel.df[keep_rows,]
head(rel_VIU_F2.df)

# Keep only those that have VIU_F2 in both
rel_VIU_F2.df <- separate(data = rel_VIU_F2.df, col = "ind1.id", into = c("pop1", "ind1"), sep = "__", remove = F)
rel_VIU_F2.df <- separate(data = rel_VIU_F2.df, col = "ind2.id", into = c("pop2", "ind2"), sep = "__", remove = F)
head(rel_VIU_F2.df)
rel_VIU_F2.df$selection_factor  <- paste0(rel_VIU_F2.df$pop1, "__", rel_VIU_F2.df$pop2)
head(rel_VIU_F2.df)
rel_VIU_F2.df <- rel_VIU_F2.df[rel_VIU_F2.df$selection_factor=="VIU_F2__VIU_F2", ]

head(rel_VIU_F2.df)
rel_VIU_F2.df <- separate(data = rel_VIU_F2.df, col = "ind1", into = c("family1", "individual1"), sep = "-", remove = F)
rel_VIU_F2.df <- separate(data = rel_VIU_F2.df, col = "ind2", into = c("family2", "individual2"), sep = "-", remove = F)
head(rel_VIU_F2.df)

fams <- unique(rel_VIU_F2.df$family1)
slice <- NULL; result.vec <- NULL; min.vec <- NULL
for(i in 1:length(fams)){
  
  slice <- rel_VIU_F2.df[rel_VIU_F2.df$family1==fams[i] & rel_VIU_F2.df$family2==fams[i], ]
  
  print(paste0(fams[i], " mean Ritland is ", round(mean(slice$ritland), digits = 2)))
  result.vec <- c(result.vec, mean(slice$ritland))
  
  print(paste0(fams[i], " min Ritland is ", round(min(slice$ritland), digits = 2)))
  min.vec <- c(min.vec, min(slice$ritland))
  
}

mean(result.vec) # 0.431
mean(min.vec).   # 0.259
length(result.vec)

# Set cutoffs based on the above
cutoff <- 0.26
#cutoff <- 0.43

# note: need to source from dev scripts
id_close_kin(cutoff = cutoff, statistic = "ritland")
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

# Remove cultured CHN samples, previously unknown but now noted as QDC samples
QDC.inds
drop.inds <- c(drop.inds, QDC.inds)
drop.inds <- unique(drop.inds)

keep.inds <- setdiff(indNames(obj), drop.inds)

obj_purged_relatives <- obj[(keep.inds)]
table(pop(obj_purged_relatives))

# And before removing? 
keep.inds <- setdiff(indNames(obj), QDC.inds)
obj_wo_QDC <- obj[(keep.inds)]
table(pop(obj_wo_QDC))


#### 04. Analysis ####
## Multivariate
# PCA from genind
pca_from_genind(data = obj_purged_relatives, PCs_ret = 4, colour_file = "00_archive/formatted_cols.csv")

# # DAPC from genind
# dapc_from_genind(data = obj, plot_allele_loadings = TRUE, colour_file = "00_archive/formatted_cols.csv")

## Dendrogram
#make_tree(bootstrap = TRUE, boot_obj = obj_purged_relatives, nboots = 10000, dist_metric = "edwards.dist", separated = FALSE)

## Drop low sample size pops
table(pop(obj_purged_relatives))
drop_pops(df = obj_purged_relatives, drop_by_pop_size = TRUE, min_indiv = 15)
obj_pop_filt
table(pop(obj_pop_filt))

## Genetic differentiation
calculate_FST(format = "genind", dat = obj_pop_filt, separated = FALSE, bootstrap = TRUE)


##### Additional analysis of the per loc stats related to the reason for selecting markers ####
obj # should be 312 inds, 409 loci
per_locus_stats(data = obj)

head(per_loc_stats.df)
dim(per_loc_stats.df)

# Load the reason for designing each locus file
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
dim(reason.df)

# Combine
all.df <- merge(x = per_loc_stats.df, y = reason.df, by = "mname", all.x = T)
head(all.df)
dim(all.df)
multi_reason_markers.vec <- all.df[duplicated(x = all.df$mname), "mname"] # these markers were selected for multiple reasons, both top_FST and top_HOBS

all.df[all.df$mname %in% multi_reason_markers.vec, ]
all.df[all.df$mname %in% multi_reason_markers.vec, "reason"] <- "FST_and_HOBS"

# Drop duplicated, as now covered in both
all.df <- all.df[!duplicated(x = all.df$mname), ]
dim(all.df)

table(all.df[all.df$Hobs > 0.25, "reason"])
length(all.df[all.df$Hobs > 0.25, "reason"])

table(all.df[all.df$Hobs <= 0.25, "reason"])
length(all.df[all.df$Hobs <= 0.25, "reason"])




save.image(file = "03_results/popgen_after_purging_close_relatives.RData")




