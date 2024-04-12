## Analysis script for OCP23 (parentage analysis), Part 2
# Ben Sutherland and Liam Surry, VIU
# initialized 2023-08-24

# Requires that OCP23_analysis_2023-08-28.R was previously run

# Clear space, then source simple_pop_stats
# Load the data
load(file = "03_results/post-filters_prepared_for_parentage.RData")

# PCA
pca_from_genind(data = obj, PCs_ret = 4, plot_eigen = F, plot_allele_loadings = F, retain_pca_obj = T)

#### 04. Prepare rubias input file ####
## Convert genepop to Rubias format

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
datatype <- "SNP" # required for genepop_to_rubias_SNP
# this is for annotate_rubias(), for an unknown reason it requires the name micro_stock_code.FN

## Convert genepop to rubias
as.data.frame(cbind(indNames(obj), as.character(pop(obj)))) 
obj

indNames(obj) # note: still in the OCP individual randomized identifier format

## Rename individuals to the alternate identifier
# Use annotated df
head(indiv_annot.df)
indiv_annot.df <- indiv_annot.df[,c("indiv", "alt.ID")]

# Bring out the individual names from the obj
indiv <- indNames(obj)
indiv <- as.data.frame(indiv)
head(indiv)

# Rename the existing samples with the information from the annotation file
indiv.df <- merge(x = indiv, y = indiv_annot.df, by = "indiv", all.x = T, sort = F)
cbind(indiv, indiv.df)

# Rename using the alt identifier
indNames(obj) <- indiv.df$alt.ID

# Remove the atypical character
indNames(obj) <- gsub(pattern = "/", replacement = "_", x = indNames(obj))
indNames(obj)

# Now that the individuals have been renamed, the original ind-to-pop pop map will no longer work
# so we need to create a new one as follows
renamed_pop_map.df <- cbind(indNames(obj), as.character(pop(obj)))
renamed_pop_map.df <- as.data.frame(renamed_pop_map.df)
colnames(renamed_pop_map.df) <- c("indiv", "pop") 
head(renamed_pop_map.df)

write.table(x = renamed_pop_map.df, file = "00_archive/renamed_ind-to-pop.txt", sep = "\t", quote = F
            , row.names = F, col.names = T
)

## All filtered loci: write to rubias for parentage assignment
pop_map.FN <- "00_archive/renamed_ind-to-pop.txt" # renamed samples
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

# Go to OCP23_analysis_part_3_2024-02-23.R
