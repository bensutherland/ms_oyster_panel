# Script to analyze families for number polymorphic markers per parental pair
# B. Sutherland, initialized 2023-12-21

# Clear space, then source simple_pop_stats start script


#### 01. Load previous data or other metadata ####
load("03_results/completed_popgen_analysis.RData")


#### 02. Polymorphic loci per parental pair ####
# Load pedigree family map
family_map.df <- read.table(file = "02_input_data/family_map.csv", header = T, sep = ",")
family_map.df

# Drop families that do not have both parents
family_map.df <- family_map.df[family_map.df$parent.1 %in% indNames(obj),]
family_map.df <- family_map.df[family_map.df$parent.2 %in% indNames(obj),]

family_map.df

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
dim(num_fams_polymorph.df)
