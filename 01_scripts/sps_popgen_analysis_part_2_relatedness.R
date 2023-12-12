# Relatedness analysis script and post-purged sibs analyses for the pilot study
# Sutherland Bioinformatics
# Initialized 2022-09-18

#### 00. Set up ####
# Source simple_pop_stats and choose Pacific oyster
load(file = "03_results/output_post-filters.Rdata")

# Working with the following: 
obj

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

# Remove cultured CHN samples, previously unknown but now noted as QDC samples
QDC.inds
drop.inds <- c(drop.inds, QDC.inds)
drop.inds <- unique(drop.inds)

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
table(pop(obj_pop_filt))

## Genetic differentiation
calculate_FST(format = "genind", dat = obj_pop_filt, separated = FALSE, bootstrap = TRUE)
