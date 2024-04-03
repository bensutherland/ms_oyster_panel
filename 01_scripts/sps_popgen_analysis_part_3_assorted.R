# Additional analysis script and post-purged sibs analyses for the pilot study
# Sutherland Bioinformatics
# Initialized 2022-09-18

#### 00. Set up ####
# Source simple_pop_stats and choose Pacific oyster
load(file = "03_results/output_post-filters.Rdata")

# Working with the following: 
obj

####### Polymorphic per pop #####
# Separate genind into individual pop genind
my_pops.list <- seppop(x = obj, drop=TRUE)
maf_cutoff <- 0.1

pop_of_interest <- NULL; poly_loci.list <- list(); poly_loci_maf.list <- list()
for(i in 1:length(my_pops.list)){
  
  # Select the population, set variable
  pop_of_interest <- names(my_pops.list)[i]
  
  # Drop monomorphic loci and store nloc in list
  print("Dropping non-polymorphic loci")
  drop_loci(df = my_pops.list[[i]], drop_monomorphic = TRUE)
  poly_loci.list[[pop_of_interest]] <- nLoc(obj_filt)
  
  # Drop under MAF threshold and store nloc in list
  print(paste0("Dropping loci below MAF ", maf_cutoff))
  maf_filt(data = my_pops.list[[i]], maf = maf_cutoff)
  
  poly_loci_maf.list[[pop_of_interest]] <- nLoc(obj_maf_filt)
  
  
}

# View results
poly_loci.list
poly_loci_maf.list


