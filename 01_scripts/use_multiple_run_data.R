# Analyze multiple runs together, where the runs may contain the same individuals
# Check consistency as well as determine the best technical replicate based on number of amplicons typed
# Sutherland Bioinformatics
# Initialized 2022-11-09

# Source simple_pop_stats and choose Pacific oyster
options(scipen = 999999999)

#### 00. Load Data ####
# Using the prompt below, load both of the following, and then save as individually-named objects
#"../amplitools/03_results/R_2022_08_04_S5XL.xls_genetic_data_only_final.gen"
#"../amplitools/03_results/R_2022_10_07_S5XL.xls_genetic_data_only_final.gen"

load_genepop(datatype = "SNP") 

# obj.r_2022_08_04 <- obj
# obj.r_2022_10_07 <- obj

# What markers are common to both obj? This is required for repooling
common_loci <- intersect(x = locNames(obj.r_2022_08_04), y = locNames(obj.r_2022_10_07))

nLoc(obj.r_2022_08_04)
obj.r_2022_08_04.common <- obj.r_2022_08_04[loc=common_loci]
nLoc(obj.r_2022_08_04.common)

nLoc(obj.r_2022_10_07)
obj.r_2022_10_07.common <- obj.r_2022_10_07[loc=common_loci]
nLoc(obj.r_2022_10_07.common)

obj <- repool(obj.r_2022_08_04.common, obj.r_2022_10_07.common)

indNames(obj) # note: Run 2022_10_07 has technical replicates in the 600 and 900 series barcodes (two copies of all samples). These are also present on 2022_08_04 run. 

#### 01. Find the best individual from the multiple runs ####
# Percent missing per individual
percent_missing_by_ind(df = obj)

#head(missing_data.df)
rownames(missing_data.df) <- seq(1:nrow(missing_data.df))
head(missing_data.df)
dim(missing_data.df)

# Separate identifier column into the component parts
missing_data.df <- separate(data = missing_data.df, col = "ind", into = c("run", "barcode", "indiv")
                            , remove = F, sep = "__")
head(missing_data.df)


# Per individual, find the one with the highest typed markers
indiv <- unique(missing_data.df$indiv)
length(indiv) # 370 unique indivs

ioi <- NULL; slice <- NULL; keep.vec <- NULL; keep <- NULL
for(i in 1:length(indiv)){
  
  # select indiv
  ioi <- indiv[i]
  slice <- missing_data.df[missing_data.df$indiv==ioi, ]
  
  # Put in descending order of ind.num.typed
  slice <- slice[order(slice$ind.num.typed, decreasing = TRUE), ]
  
  # Identify the best ind to keep (the top row)
  keep <- slice[1,"ind"]
  
  keep.vec <- c(keep.vec, keep)
  
}

keep.vec

# Retain only the best from the obj
obj.best <- obj[i=keep.vec]

# Clear the population attribute
pop(obj.best) <- rep("unkn", times = nInd(obj.best))
obj.best

# These are the ind names
indNames(obj.best)
# obj.best can now be used as an input to "simple_pop_stats_Cgig_analysis_2022-09-12.R"



#### 02. Compare technical replicates within the run ####
# Here we will use only the 2022_10_07 data, as this has the best technical replicates included
obj.r_2022_10_07

# Convert to a df
obj.df <- genind2df(obj.r_2022_10_07)

obj.df[1:5, 1:5]
colnames(obj.df)
obj.df$indiv <- rownames(obj.df)
obj.df[1:5, 577:587]

# Separate the indiv ID into components
obj.df <- separate(data = obj.df, col = "indiv", into = c("run", "barcode", "indiv"), sep = "__", remove = F)
table(obj.df$indiv) # note that all samples have two entries, this is the 600-level and 900-level barcode

# Identify the indiv IDs that have a technical replicate present
tech_rep_indivs <- dimnames(table(obj.df$indiv)==2)[[1]] 


# ### Testing out the theory for the loop ###
# test <- obj.df[obj.df$indiv=="BR1.12-A", ]
# dim(test)
# test[1:2, 1:10]
# 
# table(test[1,]==test[2,]) # Note that this ignores any instances of NAs
# 
# # Match
# test[,"96509"]
# test[1,"96509"]==test[2,"96509"]
# 
# # Different genotype (het vs. homo alt)
# test[,"714071"]
# test[1,"714071"]==test[2,"714071"]
# 
# # No data
# test[,"395635"]
# test[1,"395635"]==test[2,"395635"] # If both entries are NA, the result will be NA
# 
# # Data for only one replicate
# test[,"100388"]
# test[1,"100388"]==test[2,"100388"] # If one of the entries are NA, the result is an NA
# 
# ### /END/ Testing out the theory for the loop ###


# Loop to count matching genotypes per marker per individual
soi <- NULL; slice <- NULL; result.list <- list(); all_result.df <- NULL; 
num_false <- NULL; num_true <- NULL; percent_true <- NULL; num_typed_in_both <- NULL

for(i in 1:length(tech_rep_indivs)){
  
  print(i)
  
  # Identify the sample of interest
  soi <- tech_rep_indivs[i]
  
  # Obtain the data for the sample of interest
  slice <- obj.df[obj.df$indiv==soi, ]
  slice <- slice[grep(pattern = "run|barcode|indiv|pop", x = colnames(slice), invert = T)] # exclude those cols

  # Does the first sample geno equal the second sample geno?
  num_false <- table(slice[1,]==slice[2,])[1]
  num_true <- table(slice[1,]==slice[2,])[2]
  
  # What is the percentage? 
  num_typed_in_both <- (num_true + num_false)
  percent_match <- num_true / num_typed_in_both
  
  # Make the result into a dataframe and give the column name as the sample of interest
  result.df <- as.data.frame(c(num_false, num_true, percent_match, num_typed_in_both))
  colnames(result.df) <- soi
  
  rownames(result.df) <- c("num_false", "num_true", "percent_match", "num_typed_in_both")
  
  # Build the output df
  if(i==1){
    
    all_result.df <- result.df
    
  }else if(i > 1){
    
    all_result.df <- cbind(all_result.df, result.df)
      
  }
}


all_result.df <- t(all_result.df)
head(all_result.df)
all_result.df <- as.data.frame(all_result.df)
str(all_result.df)

pdf(file = "03_results/matched_genos.pdf", width = 7, height = 5)
plot(x = all_result.df$num_typed_in_both, y = all_result.df$percent_match, ylim = c(0,1)
     , ylab = "Percentage of matching genotypes"
     , xlab = "Number loci typed in both technical replicates")
# summary(all_result.df$percent_match)["Mean"]
# mean(all_result.df$percent_match, na.rm = T)

med.match <- median(all_result.df$percent_match, na.rm = T)
text(paste0("median = ", round(med.match, digits = 3)), x = 400, y = 0.4)
dev.off()

write_delim(x = all_result.df, file = "03_results/matched_genos.txt", delim = "\t")
