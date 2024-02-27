# Script to analyze families to inspect for null alleles
# B. Sutherland, initialized 2023-12-21

# Requires that input scripts have been run 
# Clear space and load sps
load(file = "03_results/post-filters_prepared_for_parentage_rubias_built.RData")

# Set variable depending on the dataset
#report.FN <- "../amplitools_OCP23_v.0.3/03_results/ckmr_input_rubias_135_ind_343_loc_2024-02-26_F1_vs_F0_2024-02-26/po_F0_vs_F1_pw_logl_5_report.txt"
report.FN <-  "../amplitools_v.0.5_2024-02-20/03_results/ckmr_input_rubias_142_ind_380_loc_2024-02-27_VIU_F2_vs_VIU_F1_2024-02-27/po_VIU_F1_vs_VIU_F2_pw_logl_5_report.txt"  # pilot study

#pop_map.FN <- "00_archive/renamed_ind-to-pop.txt" # OCP study
pop_map.FN <- "02_input_data/my_data_ind-to-pop_annot.txt" # pilot study

# Note: make sure you have already set the following variables
report.FN


#### 00. Identify strong empirical trios as 'true' relationships ####
### ID empirically-identified trios
# Load data
po_report.df <- read.table(file = report.FN, header = T, sep = "\t")
head(po_report.df)
po_report.df$p1_logl <- as.numeric(po_report.df$p1_logl)
po_report.df$p2_logl <- as.numeric(po_report.df$p2_logl)

## Limit to only those strong IDs
# Remove missing data
po_report.df <- po_report.df[!is.na(po_report.df$p1_logl) & !is.na(po_report.df$p2_logl), ] 

# logl >= 10 in both parents
po_report.df <- po_report.df[po_report.df$p1_logl >= 10 & po_report.df$p2_logl >= 10 & po_report.df$other_assigns=="", ]
dim(po_report.df)
head(po_report.df)

## Define empirical family map
# Who are the unique parents from these high quality assignments, 
#   sort them and create an empirical sorted family ID
parents <- NULL; parent.vec <- NULL
for(i in 1:nrow(po_report.df)){
  
  parents <- sort(c(po_report.df$p1[i], po_report.df$p2[i]))
  
  parents <- paste0(parents, collapse = "__")
  
  # Build the parent vector
  parent.vec <- c(parent.vec, parents)
  
}

empirical_family_map.df <- as.data.frame(parent.vec)
colnames(empirical_family_map.df)[1] <- "family.id"
head(empirical_family_map.df)

# Add other cols
empirical_family_map.df$offspring <- po_report.df$indiv
empirical_family_map.df <- separate(data = empirical_family_map.df, col = "family.id", into = c("parent1", "parent2"), sep = "__", remove = F)

empirical_family_map.df <- empirical_family_map.df[,c("offspring", "family.id", "parent1", "parent2")]
head(empirical_family_map.df)

# How many times were each set of parents observed together? 
family_freq.df <- as.data.frame(table(empirical_family_map.df$family.id))
# Only keep those with more than 1 observation
family_freq.df <- family_freq.df[family_freq.df$Freq > 1, ]
family_freq.df
select_empirical_families <- family_freq.df$Var1

# Keep only those families with more than one occurrences, this will comprise the empirical trios
empirical_family_map.df <- empirical_family_map.df[empirical_family_map.df$family.id %in% select_empirical_families, ]
empirical_family_map.df
dim(empirical_family_map.df)
length(unique(empirical_family_map.df$family.id))
unique(empirical_family_map.df$family.id)


#### 01.a. Back up prepared obj ####
obj.bck <- obj


#### 01.b. Retain only the empirically-defined trios in the dataset ####
offspring_to_keep <- empirical_family_map.df$offspring
parents_to_keep   <- unique(c(empirical_family_map.df$parent1,  empirical_family_map.df$parent2))
inds_to_keep <- c(offspring_to_keep, parents_to_keep)
length(offspring_to_keep) # 36 (pilot), 41 (OCP)
length(parents_to_keep)   # 15 (pilot), 15 (OCP)
length(inds_to_keep)      # 51 (pilot), 56 (OCP)

obj <- obj[inds_to_keep]


#### 02. Determine expected offspring genos ####

# What are the unique parent combinations? 
empirical_family_map.df
unique_families.df <- empirical_family_map.df[!duplicated(x = empirical_family_map.df$family.id), ]
unique_families.df <- unique_families.df[, c("family.id", "parent1", "parent2")]

# Characterize the potential offspring genos per family
# Set nulls
p1 <- NULL; p2 <- NULL; fam <- NULL; family_marker.list <- list()
for(i in 1:nrow(unique_families.df)){
  
  # Set variables
  fam <- unique_families.df$family.id[i]
  fam <- as.character(fam) # necessary in case it is being interpreted as an integer
  p1 <-  unique_families.df$parent1[i]
  p2 <-  unique_families.df$parent2[i]
  
  # Reporting
  print(paste0("Analyzing empirical family ", fam, ", comprised of parents: ", p1, " and ", p2))
  
  # Subset to only the parents for the family
  parents.obj <- obj[c(p1, p2)]
  
  # Reporting
  print("Determining offspring expected genos")
  
  # View parent genotypes, two-column data (i.e., each marker has .01 and .02, and a count of the number of alleles for each)
  parents.obj$tab[,1:5]
  
  # Drop the second call
  single_genos.df <- parents.obj$tab[, grep(pattern = "\\.02", x = colnames(x = parents.obj$tab), invert = T)]
  single_genos.df[,1:10] # only .01 is now present
  
  # Edit geno matrix per parent, converting numeric to genotype calls, derived from the first column call
  for(n in 1:nrow(single_genos.df)){
    
    single_genos.df[n,] <- gsub(pattern = "1", replacement = "het", x = single_genos.df[n,])
    single_genos.df[n,] <- gsub(pattern = "2", replacement = "homo.ref", x = single_genos.df[n,])
    single_genos.df[n,] <- gsub(pattern = "0", replacement = "homo.alt", x = single_genos.df[n,])
    
  }
  
  single_genos.df[,1:10]
  
  # Note: as shown below, there are missing values present, so need to flag these out later in the analysis
  table(single_genos.df, useNA = "ifany")
  
  # Observe the new format
  single_genos.df[,1:10]
  
  # Save out parental genos per marker
  write.csv(x = single_genos.df, file = paste0("03_results/parental_genos_per_marker_", fam, ".csv"), quote = F
              )
  
  ## Determine expected offspring genos
  # Produce matrix of expected genos per marker in offspring
  mname <- NULL; opts <- NULL; marker.list <- list()
  for(m in 1:ncol(single_genos.df)){
    
    # Debugging
    #print(m)
    
    # Identify marker name for this round
    mname <- colnames(single_genos.df)[m]
    mname
    
    # Extract the parental genotypes for this marker
    opts <- single_genos.df[,mname]
    opts
    
    # If either of the parents are NA, then need to flag the row
    if(sum(is.na(opts)) > 0){
      
      # Missing data, flag it as such
      opts <- "missing"
      
    }else if(sum(is.na(opts)) == 0){
      
      # No missing data, continue
      
    }
    
    # Simplify the parental genotypes by sorting and combining
    opts <- sort(unique(opts)) # note: order does not matter
    opts <- paste0(opts, collapse = "_")
    opts
    
    # Based on the simplified parental genotypes, what are the possible offspring genotypes? 
    if(opts=="missing"){
      
      # At least one parent is missing, so flag expected offspring as such
      outcome <- "missing"
      
    }else if(opts=="het_homo.ref"){
      
      # Parents are Aa x AA
      outcome <- c("het", "homo.ref")
      
    }else if(opts=="homo.alt"){
      
      # Parents are aa x aa
      outcome <- c("homo.alt")
      
    }else if(opts=="het"){
      
      # Parents are Aa x Aa
      outcome <- c("homo.ref", "het", "homo.alt")
      
    }else if(opts=="homo.alt_homo.ref"){
      
      # Parents are aa x AA
      outcome <- c("het")
      
    }else if(opts=="homo.ref"){
      
      # Parents are AA x AA
      outcome <- c("homo.ref")
      
    }else if(opts=="het_homo.alt"){
      
      # Parents are Aa x aa
      outcome <- c("het", "homo.alt")
      
    }else{
      
      stop("Error, outcome not designated")
      
    }
    
    marker.list[[mname]] <- outcome
    
    marker.list
  }
  
  print("Adding family to the family_marker.list")
  family_marker.list[[fam]] <- marker.list
  
}

names(family_marker.list)
length(family_marker.list)
#table(unlist(family_marker.list[[1]])) # Note NAs are still present

# Confirm all loci are present 
for(i in 1:length(family_marker.list)){
  
  print(length(family_marker.list[[i]]))
  
}

# The output is a list, where each slot is the expected offspring genotypes per locus based on the parental genotypes, 
#  and the number of slots is based on the number of families
head(family_marker.list[[1]])


#### 03. Evaluate offspring against parents ####
# Who are the offspring of the family? 
indNames(obj) # individuals available
empirical_family_map.df

families <- unique(empirical_family_map.df$family.id)

fam <- NULL; p1 <- NULL; p2 <- NULL; keep <- NULL; offspring <- NULL
single_genos_eval.df <- NULL
for(o in 1:length(families)){
  
  # Set family
  fam <- families[o]
  p1 <- gsub(pattern = "__.*", replacement = "", x = fam)
  p2 <- gsub(pattern = ".*__", replacement = "", x = fam)
  #off.string <- family_map.df[o,"offspring.string"]
  
  # Reporting
  print(paste0("Analyzing family ", fam, ", comprised of parents: ", p1, " and ", p2))
  print(paste0("Empirically defined offspring are selected from the empirical family map"))
  
  # Retain only the offspring empirically determined to be offspring of this specific family
  keep <- empirical_family_map.df[empirical_family_map.df$family.id==fam, "offspring"]
  offspring <- obj[keep]
  print(indNames(offspring))
  
  #head(family_marker.list[[fam]], n = 10)
  #offspring$tab[,1:10]
  
  # Create empty df
  single_genos.df <- NULL
  
  # For the offspring genotypes, remove the second allele column, keeping only the first allele column
  single_genos.df <- offspring$tab[, grep(pattern = "\\.02", x = colnames(x = offspring$tab), invert = T)]
  #colnames(single_genos.df)
  
  # For these offspring, edit the single_genos df to give the genotype as calculated from the single col data
  for(n in 1:nrow(single_genos.df)){
    
    single_genos.df[n,] <- gsub(pattern = "1", replacement = "het", x = single_genos.df[n,])
    single_genos.df[n,] <- gsub(pattern = "2", replacement = "homo.ref", x = single_genos.df[n,])
    single_genos.df[n,] <- gsub(pattern = "0", replacement = "homo.alt", x = single_genos.df[n,])
    
  }
  
  #single_genos.df[,1:10]
  
  # Write out offspring geno calls
  write.csv(x = single_genos.df, file = paste0("03_results/offspring_geno_calls_", fam, ".csv"))
  
  ## Check the observed offspring genotypes against the expected offspring genotypes, derived from the parents
  # Prepare a T/F matrix to fill
  single_genos_eval.df <- single_genos.df
  
  mname <- NULL
  for(i in 1:ncol(single_genos.df)){
    
    # Check each marker
    mname <- colnames(single_genos.df)[i]
    
    # What are the accepted (expected) options for this marker?
    accepted.opts <- unlist(family_marker.list[[fam]][mname])
    
    # Also allow NA to be an accepted outcome (i.e., it does not indicate an error)
    accepted.opts <- c(accepted.opts, NA) # allow NA in offspring to be accepted
    
    # If the parental genotype was missing, and therefore not able to be evaluated, denote it as such
    if(sum(accepted.opts %in% "missing") > 0){
      
      single_genos_eval.df[,i] <- "parent.missing"
      
    # If the parental genotypes were present, and the offspring genotypes predicted, continue
    }else if(sum(accepted.opts %in% "missing") == 0){
      
      # Are the observed genotypes expected genotypes? 
      single_genos_eval.df[,i] <- single_genos.df[,i] %in% accepted.opts
      
    }
    
    # Then add back offspring NAs as missing
    single_genos_eval.df[is.na(single_genos.df[,i]), i] <- "offspring.NA"
    
  }
  
  # # Confirm the edit for the NAs by checking before and after their removal/addition
  # table(single_genos.df, useNA = "ifany")
  # table(single_genos_eval.df, useNA = "ifany")
  
  write.csv(x = single_genos_eval.df, file = paste0("03_results/offspring_true_and_false_matches_", fam, ".csv"))
  
  
}


#### 04. Prepare final output report ####
files.FN <- list.files(path = "03_results/", pattern = "offspring_true_and_false_matches_")
files.FN <- files.FN[grep(pattern = "offspring_true_and_false_matches_all.csv", x = files.FN, invert = T)] # ignore the final output
files.FN

# Read in all of the T/F match matrices and combine them into one big dataframe
temp <- NULL; all_data.df <- NULL
for(i in 1:length(files.FN)){
  
  temp <- read.csv(file = paste0("03_results/", files.FN[i]))
  all_data.df <- rbind(all_data.df, temp)
  
}

dim(all_data.df)
all_data.df[1:5,1:10]

# Clean up names a bit
colnames(all_data.df)[which(colnames(all_data.df)=="X")] <- "indiv"
colnames(all_data.df) <- gsub(pattern = "^X", replacement = "", x = colnames(all_data.df))
all_data.df[1:5,1:40]
#str(all_data.df)

# Improve visualization
# Convert all to character
library("dplyr")
all_data.df <- all_data.df %>%
                   mutate_all(as.character)
all_data.df[1:5,1:10]

# Replace TRUE with dashes to make visually easier to view
library("tidyverse")
all_data.df <- all_data.df %>%
           mutate_all(funs(str_replace(., "TRUE", "-")))
all_data.df[1:5,1:10]

# Remove the first allele indicator
colnames(all_data.df) <- gsub(pattern = "\\.01", replacement = "", x = colnames(all_data.df))
all_data.df[1:5,1:10]


## Summarize the results in a dataframe
# Create matrix with the number of rows being the number of columns (i.e., number of loci evaluated), and four cols
# Note: the first column is allowed as a dummy
result.df <- matrix(data = NA, nrow = ncol(all_data.df), ncol = 4)
dim(result.df) 

# Set column names
colnames(result.df) <- c("mname", "exp.offsp.geno", "unexp.offsp.geno", "percent.exp.offsp.geno")
result.df <- as.data.frame(result.df) # make df

# Fill dataframe
exp.offsp.geno <- NULL; unexp.offsp.geno <- NULL; percent.exp.offsp.geno <- NULL; mname <- NULL
for(i in 1:ncol(all_data.df)){
  
  mname                 <- colnames(all_data.df)[i]
  result.df[i, "mname"] <- mname
  
  exp.offsp.geno        <- sum(all_data.df[,i]=="-", na.rm = T)
  result.df[i, "exp.offsp.geno"] <- exp.offsp.geno
  
  unexp.offsp.geno      <- sum(all_data.df[,i]=="FALSE", na.rm = T)
  result.df[i, "unexp.offsp.geno"] <- unexp.offsp.geno
  
  percent.exp.offsp.geno      <- exp.offsp.geno / (exp.offsp.geno + unexp.offsp.geno)
  result.df[i, "percent.exp.offsp.geno"] <- percent.exp.offsp.geno
  
}

# Drop the unneeded column
result.df <- result.df[result.df$mname!="indiv", ]

# Adjust decimal
result.df$percent.exp.offsp.geno <- formatC(x = result.df$percent.exp.offsp.geno, digits = 3, format = "f")

head(result.df)

# What is the distribution of the erroneous calls?
pdf(file = "03_results/hist_per_locus_num_unexpected_genos.pdf", width = 6.5, height = 4)
hist(result.df$unexp.offsp.geno
     , main = ""
     , las = 1
     , xlab = "Per locus, number of indiv. with unexpected genos"
     )
dev.off()


# How many loci have at least two erroneous calls? 
# table(result.df$unexp.offsp.geno >= 2) # 97 (pilot), 80 (OCP)
table(result.df$unexp.offsp.geno >= 4) #  (pilot), 43 (OCP)

# Write out per locus info
write.csv(x = all_data.df, file = paste0("03_results/per_locus_all_results.csv"), row.names = F)
write.csv(x = result.df, file = paste0("03_results/per_locus_expected_offsp_genos.csv"), row.names = F)


# What loci are erroneous in at least two offspring? 
# loci_to_drop <- result.df[result.df$unexp.offsp.geno >= 2, "mname"]
loci_to_drop <- result.df[result.df$unexp.offsp.geno >= 4, "mname"]

# write it out
write.table(x = loci_to_drop, file = "03_results/markers_with_unexpected_genos_in_offspring.txt"
            , quote = F, sep = "\t", row.names = F, col.names = F
            )


#### Remove the problematic loci from the original genind file #####
obj <- obj.bck
obj

drop_loci(df = obj, drop_file = "03_results/markers_with_unexpected_genos_in_offspring.txt")
obj_filt

#### Write the new data to a rubias file ####
# Write out to rubias
pop_map.FN  # renamed samples

genepop_to_rubias_SNP(data = obj_filt, sample_type = "reference", custom_format = TRUE, micro_stock_code.FN = micro_stock_code.FN
                      , pop_map.FN = pop_map.FN)

print("Your output is available as '03_results/rubias_output_SNP.txt")


#### Rename rubias file ####
# Obtain some variables to create filename
date <- format(Sys.time(), "%Y-%m-%d")
indiv.n <- nInd(obj_filt)
loci.n  <- nLoc(obj_filt)
rubias_custom.FN <- paste0("03_results/rubias_", indiv.n, "_ind_", loci.n, "_loc_", date, ".txt")

# Copy to rename rubias output file
file.copy(from = "03_results/rubias_output_SNP.txt", to = rubias_custom.FN, overwrite = T)
print(paste0("The output is saved as ", rubias_custom.FN))

# Save image
save.image("03_results/post-filters_prepared_for_parentage_rubias_built_problem_genos_rem_rubias_built.RData")

print("Here you need to copy the above rubias file to amplitools results folder.")

# Go to OCP23_analysis_part_5_2024-02-26.R or 
