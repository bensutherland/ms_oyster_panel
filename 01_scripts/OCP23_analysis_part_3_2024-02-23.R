## Analysis script for OCP23 (parentage analysis), Part 3
# Ben Sutherland and Liam Surry, VIU
# initialized 2023-08-24

# Requires that OCP23_analysis_part_2_2024-02-23.R was previously run
# Also requires that you have copied the prepared rubias file into amplitools


#### 05. All filtered loci parentage analysis ####
# Clear space, and launch amplitools initiator (i.e., 01_scripts/00_initiator.R)

# Set user variables
input_rubias.FN <- "03_results/rubias_135_ind_343_loc_2024-02-26.txt"
parent_pop <- "F0"
offspring_pop <- "F1"
cutoff <- 5

# Prepare run folder
date <- format(Sys.time(), "%Y-%m-%d")
input_rubias_short.FN <- gsub(pattern = "03_results/", replacement = "", x = input_rubias.FN)
input_rubias_short.FN <- gsub(pattern = ".txt", replacement = "", x = input_rubias_short.FN)

run_folder.FN <- paste0("03_results/ckmr_input_", input_rubias_short.FN, "_"
                        , offspring_pop, "_vs_", parent_pop,"_", date
                        )
print("Making new result folder...")
print(run_folder.FN)
dir.create(run_folder.FN)


# Run ckmr on the input file
ckmr_from_rubias(input.FN = input_rubias.FN
                 , parent_pop = "F0"
                 , offspring_pop = "F1"
                 , cutoff = 5
                 , output.dir = run_folder.FN
)

# # Filtered loci and pilot study null allele removed
# ckmr_from_rubias(input.FN = "03_results/rubias_output_SNP_filtered_and_null_pilot_drop.txt", parent_pop = "F0"
#                  , offspring_pop = "F1", cutoff = 5
# )

# Plot the output results
graph_relatives(input.FN = "03_results/ckmr_input_rubias_135_ind_343_loc_2024-02-26_F1_vs_F0_2024-02-26/po_F0_vs_F1_pw_logl_5.txt", logl_cutoff = 5
                , drop_string = "", directed = F, plot_width = 8, plot_height = 8
)

# # Before going to the next script, set some variables
# report.FN <- "03_results/parentage_with_all_filtered_loci/po_F0_vs_F1_pw_logl_5_report.txt"

save.image(paste0(run_folder.FN, "/ckmr_completed.RData"))

# Now go to 01_scripts/exploring_families.R

