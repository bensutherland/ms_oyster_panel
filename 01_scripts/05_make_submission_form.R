# Script to prepare submission form for amplicon panel
# B. Sutherland (Sutherland Bioinformatics)
# 2022-01-26

#### Front Matter ####
# Clean space
# rm(list=ls())

## Install and load packages
# install.packages("devtools")
# install.packages("hierfstat")
# install.packages("adegenet")
# install.packages("pegas")
# install.packages("dartR")

# library("devtools")
# library("hierfstat")
# library("adegenet")
# library("pegas")
# library("dartR")

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
output.dir <- "05_submission_form//"
input.dir <- "02_input_data//"
species <- "Cgig"

#### 1. Import data ####
#print(paste0("Loading data from ", input.FN))






#### Next Steps: #####
