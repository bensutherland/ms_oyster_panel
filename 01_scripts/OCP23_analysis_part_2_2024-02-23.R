## Analysis script for OCP23 (parentage analysis), Part 2
# Ben Sutherland and Liam Surry, VIU
# initialized 2023-08-24

# Requires that OCP23_analysis_2023-08-28.R was previously run

# Clear space, then source simple_pop_stats
# Load the data
load(file = "03_results/filtered_genind_before_ckmr.RData")

### TODO: HERE SHOULD DO A PCA ###
pca()

#### 04. Parentage analysis ####
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

# Copy to retain
file.copy(from = "03_results/rubias_output_SNP.txt", to = "../amplitools_OCP23_v.0.3/03_results/rubias_output_SNP_all_filtered_loci.txt", overwrite = T)

# ## Filtered loci but also with those removed due to null allele occurrences in pilot study
# drop_loci.FN <- "02_input_data/loci_to_remove_from_pilot.txt"
# drop_loci(df = obj, drop_file = drop_loci.FN)
# 
# genepop_to_rubias_SNP(data = obj_filt, sample_type = "reference"
#                       , custom_format = TRUE, micro_stock_code.FN = micro_stock_code.FN
#                       , pop_map.FN = pop_map.FN
# )
# print("Your output is available as '03_results/rubias_output_SNP.txt")
# # Copy to retain
# file.copy(from = "03_results/rubias_output_SNP.txt", to = "../amplitools/03_results/rubias_output_SNP_filtered_and_null_pilot_drop.txt", overwrite = T)

# Save out image
save.image("03_results/completed_popgen_analysis.RData")

### Next, go to OCP23_analysis_part_3_2024-02-23.R 


# ##### Update sample IDs in the rubias file ####
# # Update sample IDs in the rubias file
# rubias.df <- read.delim2(file = "03_results/rubias_output_SNP.txt", sep = "\t")
# dim(rubias.df)
# rubias.df[1:5, 1:10]
# rubias.df$repunit <- rubias.df$collection # hacky fix to whatever caused the error for the repunit being listed as the alt.id
# rubias.df[1:5, 1:10]
# 
# ## Load annotated df that was manually made earlier
# indiv_annot.df <- read.table(file = "00_archive/my_data_ind-to-pop_annot.txt"
#                              , header = T, sep = "\t"
# )
# head(indiv_annot.df)
# indiv_annot.df <- indiv_annot.df[,c("indiv", "alt.ID")]
# 
# rubias.df[1:5, 1:10]
# 
# rubias.df <- merge(x = rubias.df, y = indiv_annot.df, by = "indiv", all.x = TRUE, sort = F)
# dim(rubias.df)
# rubias.df[1:5, 699:709]
# 
# rubias.df <- rubias.df[, !colnames(rubias.df) %in% "indiv"]
# rubias.df[1:5, 1:10]
# 
# rubias.df <-  rubias.df %>% 
#                   select(alt.ID, everything())
# rubias.df[1:5, 1:10]
# colnames(rubias.df)[which(colnames(rubias.df)=="alt.ID")] <- "indiv"
# rubias.df[1:5, 1:10]
# 
# write.table(x = rubias.df, file = "03_results/rubias_output_SNP_renamed.txt", sep = "\t", row.names = FALSE)

# file.copy(from = "03_results/rubias_output_SNP.txt", to = "../amplitools/03_results/cgig_all_rubias.txt", overwrite = T)
# save.image(file = "03_results/completed_analysis_to_rubias.RData")



#### Parentage Analysis ####
# Clear space, and launch amplitools initiator (i.e., 01_scripts/00_initiator.R)

# Using this output, move to "amplitools/01_scripts/ckmr_from_rubias.R"
# All filtered loci
ckmr_from_rubias(input.FN = "03_results/rubias_output_SNP_all_filtered_loci.txt", parent_pop = "F0"
                 , offspring_pop = "F1", cutoff = 5
)

# Filtered loci and pilot study null allele removed
ckmr_from_rubias(input.FN = "03_results/rubias_output_SNP_filtered_and_null_pilot_drop.txt", parent_pop = "F0"
                 , offspring_pop = "F1", cutoff = 5
)



### New items to do: 
# check number of loci per indiv from rubias file here (amplitools), retain to connect to report
# use sex attribute within the prep report

# Generate report from the output of CKMR-sim
prep_report(relationship = "PO")



# Plot the output results
graph_relatives(input.FN = "03_results/po_F0_vs_F1_pw_logl_5.txt", logl_cutoff = 5
                , drop_string = "", directed = F, plot_width = 12, plot_height = 12
)

graph_relatives(input.FN = "03_results/fs_offsp_F1_pw_logl_5.txt", logl_cutoff = 5
                , drop_string = "", directed = F, plot_width = 12, plot_height = 12
)

graph_relatives(input.FN = "03_results/fs_parent_F0_pw_logl_5.txt", logl_cutoff = 5
                , drop_string = "", directed = F, plot_width = 12, plot_height = 12
)


