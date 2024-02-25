## Analysis script for OCP23 (parentage analysis), Part 3
# Ben Sutherland and Liam Surry, VIU
# initialized 2023-08-24

# Requires that OCP23_analysis_part_2_2024-02-23.R was previously run

#### 05. All filtered loci parentage analysis ####
# Clear space, and launch amplitools initiator (i.e., 01_scripts/00_initiator.R)

# Using this output, move to "amplitools/01_scripts/ckmr_from_rubias.R"
# All filtered loci
ckmr_from_rubias(input.FN = "03_results/rubias_output_SNP_all_filtered_loci.txt", parent_pop = "F0"
                 , offspring_pop = "F1", cutoff = 5
)

# # Filtered loci and pilot study null allele removed
# ckmr_from_rubias(input.FN = "03_results/rubias_output_SNP_filtered_and_null_pilot_drop.txt", parent_pop = "F0"
#                  , offspring_pop = "F1", cutoff = 5
# )

# Plot the output results
graph_relatives(input.FN = "03_results/po_F0_vs_F1_pw_logl_5.txt", logl_cutoff = 5
                , drop_string = "", directed = F, plot_width = 12, plot_height = 12
)

# Before going to the next script, set some variables
report.FN <- "03_results/parentage_with_all_filtered_loci/po_F0_vs_F1_pw_logl_5_report.txt"

save.image("03_results/parentage_completed.RData")

# Now go to 01_scripts/exploring_families.R
