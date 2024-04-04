## Analysis script for pilot (parentage analysis), Part 5
# Ben Sutherland, VIU
# initialized 2023-08-24

# Note: the exploring families script should have been run, and the rubias output with bad loci removed should have been 
#  copied to the amplitools results folder, with the custom rubias name

#### Second parentage analysis ####
# Clear space, and launch amplitools initiator (i.e., 01_scripts/00_initiator.R)

# Set user variables
cutoff <- 5

## OCP
#input_rubias.FN <- "03_results/rubias_135_ind_263_loc_2024-02-26.txt" # OCP
#input_rubias.FN <- "03_results/rubias_135_ind_300_loc_2024-02-26.txt" # OCP
# parent_pop <- "F0"
# offspring_pop <- "F1"

## Pilot
input_rubias.FN <- "03_results/rubias_142_ind_328_loc_2024-04-03.txt" # Pilot
parent_pop <- "VIU_F1"
offspring_pop <- "VIU_F2"

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
                 , parent_pop = parent_pop
                 , offspring_pop = offspring_pop
                 , cutoff = cutoff
                 , output.dir = run_folder.FN
)


# Plot the output results
graph_relatives(input.FN = "03_results/ckmr_input_rubias_142_ind_328_loc_2024-04-03_VIU_F2_vs_VIU_F1_2024-04-03/po_VIU_F1_vs_VIU_F2_pw_logl_5.txt"
                , logl_cutoff = 5
                , drop_string = "", directed = F
                , plot_width = 8, plot_height = 8
)


### There were near-identical individual identified, let's remove
rubias.df <- read.delim2(file = "03_results/rubias_142_ind_328_loc_2024-04-03.txt")
rubias.df[1:4, 1:10]
rubias.df <- rubias.df[rubias.df$indiv!="BR1.18-E", ]
write.table(x = rubias.df, file = "03_results/rubias_141_ind_328_loc_2024-04-03.txt"
            , quote = F, sep = "\t", row.names = F)

# Prepare the output folder
input_rubias.FN <- "03_results/rubias_141_ind_328_loc_2024-04-03.txt"
parent_pop <- "VIU_F1"
offspring_pop <- "VIU_F2"

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
                 , parent_pop = parent_pop
                 , offspring_pop = offspring_pop
                 , cutoff = cutoff
                 , output.dir = run_folder.FN
)


# Plot the output results
graph_relatives(input.FN = "03_results/ckmr_input_rubias_142_ind_328_loc_2024-04-03_VIU_F2_vs_VIU_F1_2024-04-03/po_VIU_F1_vs_VIU_F2_pw_logl_5.txt"
                , logl_cutoff = 5
                , drop_string = "", directed = F
                , plot_width = 8, plot_height = 8
)





### Update: no longer processing the F1 vs. F0 analysis, as not needed ####

# #### Add the F1 vs F0 ####
# ## Pilot, first generation, all loci
# input_rubias.FN <- "03_results/rubias_142_ind_380_loc_2024-02-27.txt" # Pilot
# parent_pop <- "VIU_F0"
# offspring_pop <- "VIU_F1"
# 
# # Prepare run folder
# date <- format(Sys.time(), "%Y-%m-%d")
# input_rubias_short.FN <- gsub(pattern = "03_results/", replacement = "", x = input_rubias.FN)
# input_rubias_short.FN <- gsub(pattern = ".txt", replacement = "", x = input_rubias_short.FN)
# 
# run_folder.FN <- paste0("03_results/ckmr_input_", input_rubias_short.FN, "_"
#                         , offspring_pop, "_vs_", parent_pop,"_", date
# )
# print("Making new result folder...")
# print(run_folder.FN)
# dir.create(run_folder.FN)
# 
# 
# # Run ckmr on the input file
# ckmr_from_rubias(input.FN = input_rubias.FN
#                  , parent_pop = parent_pop
#                  , offspring_pop = offspring_pop
#                  , cutoff = cutoff
#                  , output.dir = run_folder.FN
# )
# 
# 
# # Plot the output results
# graph_relatives(input.FN = "03_results/ckmr_input_rubias_142_ind_380_loc_2024-02-27_VIU_F1_vs_VIU_F0_2024-02-27/po_VIU_F0_vs_VIU_F1_pw_logl_5.txt"
#                 , logl_cutoff = 5
#                 , drop_string = "", directed = F
#                 , plot_width = 8, plot_height = 8
# )
# 
# 
# 
# ## Pilot, first generation, CLEANED
# 
# input_rubias.FN <- "03_results/rubias_142_ind_347_loc_2024-02-27.txt" # Pilot
# parent_pop <- "VIU_F0"
# offspring_pop <- "VIU_F1"
# 
# # Prepare run folder
# date <- format(Sys.time(), "%Y-%m-%d")
# input_rubias_short.FN <- gsub(pattern = "03_results/", replacement = "", x = input_rubias.FN)
# input_rubias_short.FN <- gsub(pattern = ".txt", replacement = "", x = input_rubias_short.FN)
# 
# run_folder.FN <- paste0("03_results/ckmr_input_", input_rubias_short.FN, "_"
#                         , offspring_pop, "_vs_", parent_pop,"_", date
# )
# print("Making new result folder...")
# print(run_folder.FN)
# dir.create(run_folder.FN)
# 
# 
# # Run ckmr on the input file
# ckmr_from_rubias(input.FN = input_rubias.FN
#                  , parent_pop = parent_pop
#                  , offspring_pop = offspring_pop
#                  , cutoff = cutoff
#                  , output.dir = run_folder.FN
# )
# 
# 
# # Plot the output results
# graph_relatives(input.FN = "03_results/ckmr_input_rubias_142_ind_347_loc_2024-02-27_VIU_F1_vs_VIU_F0_2024-02-27/po_VIU_F0_vs_VIU_F1_pw_logl_5.txt"
#                 , logl_cutoff = 5
#                 , drop_string = "", directed = F
#                 , plot_width = 8, plot_height = 8
# )






# Finished
