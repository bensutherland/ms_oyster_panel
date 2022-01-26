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


#install.packages("tidyr")

# library("devtools")
# library("hierfstat")
# library("adegenet")
# library("pegas")
# library("dartR")
library("tidyr")

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
resources.dir <- "02_input_data//"
input.dir <- "04_extract_loci//"
output.dir <- "05_submission_form//"

species <- "Cgig"

#### 1. Import data ####
#print(paste0("Loading data from ", input.FN))

# Import sequence data
selected_chr_and_seq <- read.delim2(file = "04_extract_loci/selected_chr_and_seq.txt", header = F)
colnames(selected_chr_and_seq) <- c("chr_info", "seq")
head(selected_chr_and_seq, n = 1)

selected_chr_and_seq <- separate(data = selected_chr_and_seq, col = "chr_info", into = c("chr", "pos_range"), sep = ":", remove = T)
selected_chr_and_seq <- separate(data = selected_chr_and_seq, col = "pos_range", into = c("lower_range", "upper_range")
                                 , sep = "-", remove = T)
head(x = selected_chr_and_seq, n = 1)
str(selected_chr_and_seq)
selected_chr_and_seq$lower_range <- as.numeric(selected_chr_and_seq$lower_range)
selected_chr_and_seq$upper_range <- as.numeric(selected_chr_and_seq$upper_range)
str(selected_chr_and_seq)

# Import marker data
selected_marker_info <- read.delim2(file = "04_extract_loci/vcf_selection.csv", header = F, sep = ",")
colnames(x = selected_marker_info) <- c("chr", "pos", "info", "ref", "alt")
head(selected_marker_info, n = 2)

selected_marker_info <- separate(data = selected_marker_info, col = "info", into = c("mname", "SNP_pos", "align_strand"), sep = ":", remove = T)
head(selected_marker_info, n = 2)

str(selected_marker_info)

## Merge data (note: depending on only a single variant per chromosome)
length(union(x = selected_chr_and_seq$chr, selected_marker_info$chr))
nrow(selected_chr_and_seq)
nrow(selected_marker_info)

# Not going to work... need better method
# Iteratively check each row and find the correct mname

chroi <- NULL; lower.oi <- NULL; upper.oi <- NULL; chunk.df <- NULL

for(i in 1:nrow(selected_chr_and_seq)){
  
  chroi <-    selected_chr_and_seq[i, "chr"]
  lower.oi <- selected_chr_and_seq[i, "lower_range"]
  upper.oi <- selected_chr_and_seq[i, "upper_range"]
  
  print(paste0("Looking for ", chroi, " within range ", lower.oi, "-", upper.oi))
  
  # Subset matching dataframe
  chunk.df <- selected_marker_info[which(selected_marker_info$chr==chroi), ]
  
  for(n in 1:nrow(chunk.df)){
    
    if(chunk.df$pos[n] > lower.oi && chunk.df$pos[n] < upper.oi){
      
      mname.oi <- chunk.df[n,"mname"]
      selected_chr_and_seq$mname[i] <- mname.oi
      
    }
    
  }
  
  
}


length(selected_chr_and_seq$mname)


all_data <- merge(x = selected_chr_and_seq, y = selected_marker_info, by = "mname")
all_data$left_seq  <- substr(x = all_data$seq, start =   1, stop = 200) 
all_data$ref_nuc   <- substr(x = all_data$seq, start = 201, stop = 201) # has reference allele last
all_data$right_seq <- substr(x = all_data$seq, start = 202, stop = 401) # has reference allele last

all_data$strand <- rep("NA", times = nrow(all_data))
all_data$mtype <- rep("SNP", times = nrow(all_data))
all_data$priority <- rep(2, times = nrow(all_data))



# CONFIRM matching REF/ and calculated REF/ 
# TODO #

all_data$formatted_seq <- paste0(all_data$left_seq, "[", all_data$ref, "/", all_data$alt, "]", all_data$right_seq)

head(all_data, n = 3)

all_data_final <- all_data[, c("mname", "chr.x", "pos", "pos"
                               , "ref", "alt", "strand", "mtype"
                               , "priority", "formatted_seq"
                               )]

head(all_data_final)

#### Next Steps: #####
# FINALIZE quality checking above
# then write out results
# May want to converta ll to capitals
# 


