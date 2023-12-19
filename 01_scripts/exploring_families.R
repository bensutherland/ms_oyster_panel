# Script to analyze families to inspect for null alleles
# B. Sutherland, initialized 2023-12-19

# Clear space, then source simple_pop_stats start script

#load("03_results/completed_popgen_analysis.RData")

# Load family map
family_map.df <- read.table(file = "02_input_data/family_map.csv", header = T, sep = ",")
family_map.df

# Drop families that do not have both parents
family_map.df <- family_map.df[family_map.df$parent.1 %in% indNames(obj),]
family_map.df <- family_map.df[family_map.df$parent.2 %in% indNames(obj),]

family_map.df

#### 01. Polymorphic loci per parental pair ####
# How many loci are polymorphic per parental pair, and 

# Analyze for each parental pair in the map
# Set nulls
p1 <- NULL; p2 <- NULL; fam <- NULL; per_family_loci.list <- list()
for(i in 1:nrow(family_map.df)){

  # set variables
  fam <- family_map.df$family.id[i]
  p1 <-  family_map.df$parent.1[i]
  p2 <-  family_map.df$parent.2[i]

  print(paste0("***Analyzing family ", fam, ", comprised of parents: ", p1, " and ", p2, "***"))

  parents.obj <- obj[c(p1, p2)]


  # Check number polymorphic
  drop_loci(df = parents.obj, drop_monomorphic = T)

  parents.obj <- obj_filt
  #print(parents.obj)

  #parents.obj
  
  # Save the names of polymorphic loci in the parents
  per_family_loci.list[[fam]] <- locNames(parents.obj)
  
}
# Reporting
print("Names of polymorphic loci per family in per_family_loci.list")

num_fams_polymorph.df <- sort(table(unlist(per_family_loci.list)), decreasing = T)
num_fams_polymorph.df <- as.data.frame(num_fams_polymorph.df)
colnames(num_fams_polymorph.df) <- c("mname", "freq")
head(num_fams_polymorph.df)
dim(num_fams_polymorph.df)

# Which markers are missing?
monomorphs <- setdiff(x = locNames(obj), y = num_fams_polymorph.df$mname)
monomorphs <- as.data.frame(monomorphs)
colnames(monomorphs) <- "mname"
head(monomorphs)
monomorphs$freq <- "0"
head(monomorphs)

# Combine
num_fams_polymorph.df <- rbind(num_fams_polymorph.df, monomorphs)

# TODO: this can be used to join with the incompatibilities table below


#### 02. Determine expected offspring genos ####
# Set nulls
p1 <- NULL; p2 <- NULL; fam <- NULL; family_marker.list <- list()
for(i in 1:nrow(family_map.df)){
  
# i <- 2

  # Set variables
  fam <- family_map.df$family.id[i]
  p1 <-  family_map.df$parent.1[i]
  p2 <-  family_map.df$parent.2[i]
  
  # Reporting
  print(paste0("Analyzing family ", fam, ", comprised of parents: ", p1, " and ", p2))
  
  # Subset to only the parents for the family
  parents.obj <- obj[c(p1, p2)]
  
  # Reporting
  print("Determining offspring expected genos")
  
  # Parent genotypes, two-column data (i.e., each marker has .01 and .02, and a count of the number of alleles for each)
  parents.obj$tab[,1:5]
  
  # Drop the second call
  single_genos.df <- parents.obj$tab[, grep(pattern = "\\.02", x = colnames(x = parents.obj$tab), invert = T)]
  
  # Edit matrix per parent, converting numeric to genotype calls, derived from the first column call
  single_genos.df[,1:10] # only .01 is now present
  
  # Determine the parent genotype based on the single column call
  for(n in 1:nrow(single_genos.df)){
    
    single_genos.df[n,] <- gsub(pattern = "1", replacement = "het", x = single_genos.df[n,])
    single_genos.df[n,] <- gsub(pattern = "2", replacement = "homo.ref", x = single_genos.df[n,])
    single_genos.df[n,] <- gsub(pattern = "0", replacement = "homo.alt", x = single_genos.df[n,])
    
  }
  
  # Note: there are missing values present, so should flag these out
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
    opts <- sort(unique(opts)) # note: order should not matter
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
table(unlist(family_marker.list[[1]])) # Note NAs are still present

# Confirm all loci are present 
for(i in 1:length(family_marker.list)){
  
  print(length(family_marker.list[[i]]))
  
}



#### 03. Evaluate offspring against parents ####

# Who are the offspring of the family? 
indNames(obj) # individuals available
family_map.df

for(o in 1:nrow(family_map.df)){
  
  # Set family
  fam <- family_map.df[o,"family.id"]
  p1 <- family_map.df[o,"parent.1"]
  p2 <- family_map.df[o,"parent.2"]
  off.string <- family_map.df[o,"offspring.string"]
  
  # Reporting
  print(paste0("Analyzing family ", fam, ", comprised of parents: ", p1, " and ", p2))
  print(paste0("Offspring are denoted by '", off.string,"'"))
  
  # Retain only the offspring denoted by the indicated string in the family map
  keep <- indNames(obj)[grep(pattern = off.string, x = indNames(obj))]
  offspring <- obj[keep]
  print(indNames(offspring))
  
  #head(family_marker.list[[fam]], n = 10)
  #offspring$tab[,1:10]
  
  single_genos.df <- NULL
  
  # Remove the second allele column, keep only the single allele column for the offspring
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
  
  
  
  ## Check observed offspring against expected offspring genos
  # Prepare a T/F matrix to fill
  single_genos_eval.df <- single_genos.df
  
  mname <- NULL
  for(i in 1:ncol(single_genos.df)){
    
    # Check each marker at a time 
    mname <- colnames(single_genos.df)[i]
    
    # What are the accepted options for this marker? 
    accepted.opts <- unlist(family_marker.list[[fam]][mname])
    
    accepted.opts <- c(accepted.opts, NA) # allow NA in offspring to be accepted
    
    # If the parental genotype was not evaluable, denote as such
    if(sum(accepted.opts %in% "missing") > 0){
      
      single_genos_eval.df[,i] <- "parent.missing"
      
    # If the parental genotype was evaluable and the offspring genos predicted, continue
    }else if(sum(accepted.opts %in% "missing") == 0){
      
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

# Read in all of the T/F match matrices and combine them into one big df
temp <- NULL; all_data.df <- NULL
for(i in 1:length(files.FN)){
  
  temp <- read.csv(file = paste0("03_results/", files.FN[i]))
  all_data.df <- rbind(all_data.df, temp)
  
}

dim(all_data.df)
all_data.df[1:5,1:10]

# Clean up a bit
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

library("tidyverse")
all_data.df <- all_data.df %>%
           mutate_all(funs(str_replace(., "TRUE", "-")))
all_data.df[1:5,1:10]


## Add a bottom row indicating the number of polymorphic
colnames(all_data.df) <- gsub(pattern = "\\.01", replacement = "", x = colnames(all_data.df))
all_data.df[1:5,1:10]


head(num_fams_polymorph.df)

num_fams_polymorph_t.df <- t(num_fams_polymorph.df)
num_fams_polymorph_t.df[,1:5]
colnames(num_fams_polymorph_t.df) <- num_fams_polymorph_t.df["mname",]
num_fams_polymorph_t.df[,1:5]

dim(num_fams_polymorph_t.df)
num_fams_polymorph_t.df[,1:5]
dim(all_data.df)
all_data.df[1:5,1:5]

# Add indiv col
dummy.col <- as.data.frame(c("NA", "NA"))
colnames(dummy.col) <- "indiv"
rownames(dummy.col) <- c("mname", "freq")
dummy.col

num_fams_polymorph_t.df <- cbind(num_fams_polymorph_t.df, dummy.col)
num_fams_polymorph_t.df

# Combine
all_data_and_num_fam_poly.df <- rbind(all_data.df, num_fams_polymorph_t.df)
dim(all_data_and_num_fam_poly.df)
all_data_and_num_fam_poly.df[,1:5]



write.csv(x = all_data_and_num_fam_poly.df, file = paste0("03_results/offspring_true_and_false_matches_all.csv"), row.names = F)



# all_data_and_num_fam_poly_t.df <- t(all_data_and_num_fam_poly.df)
# dim(all_data_and_num_fam_poly_t.df)
# all_data_and_num_fam_poly_t.df[1:5,1:5]
# write.csv(x = all_data_and_num_fam_poly_t.df, file = paste0("03_results/offspring_true_and_false_matches_all_t.csv"), row.names = T)

# for(i in 2:ncol(all_data.df)){
#   
#   print(sum(all_data.df[,i], na.rm = T))
#   
# }


