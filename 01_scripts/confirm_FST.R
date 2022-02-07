# Confirms that elevated FST markers result in an elevated FST value for the populations
# Depends on 01_identifying_markers.R already having been ran
# 2022-02-07

#### Front Matter ####
# Clean space
# rm(list=ls())

## Install and load packages
# install.packages("devtools")
# install.packages("hierfstat")
# install.packages("adegenet")
# install.packages("pegas")
# install.packages("dartR")

library("devtools")
library("hierfstat")
library("adegenet")
library("pegas")
library("dartR")

# Set working directory
if(Sys.info()["nodename"] == "Wayne.local"){ 
  print("On Wayne, ready to go")
  setwd("~/Documents/00_sutherland_bioinformatics/GBMF_UBC_Pacific_oyster/ms_oyster_panel/") 
} else if(Sys.info()["nodename"] == "Xavier"){
  print("On Xavier, need to set path")
  #setwd("~/Documents/01_moore_oyster_project/stacks_workflow_all_data/") # Xavier
} else {
  print("You are on an unrecognized system, please set working directory manually")
}


## Info
# sessionInfo()

# Set variables
output.dir <- "03_marker_selection//"

#### 1. Find the names of the top 300 FST markers
sel_mnames.df <- read.csv(file = "04_extract_loci/selected_mnames.csv", header = F)
sel_mnames.df <- as.data.frame(sel_mnames.df, stringsAsFactors = F)
sel_mnames.df$V1 <- as.character(sel_mnames.df$V1)
sel_mnames.df$V2 <- as.character(sel_mnames.df$V2)
head(sel_mnames.df)
str(sel_mnames.df)

colnames(sel_mnames.df) <- c("mname", "mtype")
head(sel_mnames.df)

# Select only those that were selected for top_FST characteristic
sel_mnames_FST.df <- sel_mnames.df[sel_mnames.df$mtype==" top_FST", ]
dim(sel_mnames_FST.df)


#### 2. Load previous results ####
load(file = paste0(output.dir, "adegenet_output.RData"))

# Main, all loci, BC indiv
target_pops_maf_filt.gid

# Subset only FST loci, BC indiv
target_pops_maf_filt_FST_only.gid <- target_pops_maf_filt.gid[loc=sel_mnames_FST.df$mname]

# FST on each
# Convert genind to hierfstat
# All data
all.data.hf <- genind2hierfstat(target_pops_maf_filt.gid)
#rownames(target_pops_maf_filt.gid) <- indNames(target_pops_maf_filt.gid) # May not be needed

# FST only markers
FST_only_markers.hf <- genind2hierfstat(dat = target_pops_maf_filt_FST_only.gid)
dim(FST_only_markers.hf)


#### 3. Pairwise Fst #####
# All data
pairwise.wc.fst <- pairwise.WCfst(all.data.hf)
write.csv(x = pairwise.wc.fst, file = paste0(output.dir, "wcfst_all_data.csv"))

# Selected data
pairwise.wc.fst_top_FST <- pairwise.WCfst(FST_only_markers.hf)
write.csv(x = pairwise.wc.fst_top_FST, file = paste0(output.dir, "wcfst_top_FST_markers.csv"))

# Bootstrapping
nboots <- 1000

# Select which dataset to run: 
bootstrapped_input.hf <- all.data.hf
output.FN <- "wcfst_all_data_nboot_"

# OR
bootstrapped_input.hf <- FST_only_markers.hf
output.FN <- "wcfst_top_FST_nboot_"

# Calculations on the selected input
boot.fst.all <- boot.ppfst(dat = bootstrapped_input.hf, nboot = 1000, quant = c(0.025,0.975))
boot.fst.all
# note that nboot = 1000 is about the same as nboot = 10,000 (very marginally different)

# Collect output
lower.limit <- t(boot.fst.all$ll)
upper.limit <- boot.fst.all$ul
upper.limit[is.na(upper.limit)] <- 0
lower.limit[is.na(lower.limit)] <- 0
boot.fst.all.output <- upper.limit + lower.limit
boot.fst.all.output

filename <- paste0(output.dir, output.FN, nboots, ".csv")
write.csv(x = boot.fst.all.output, file = filename)

