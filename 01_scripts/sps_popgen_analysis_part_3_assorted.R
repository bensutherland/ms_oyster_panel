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


## Private alleles
regional_obj <- obj

# Combine related pops to query private alleles at regional level
unique(pop(regional_obj))
pop(regional_obj) <- gsub(pattern = "VIU_F0|VIU_F1|VIU_F2", replacement = "VIU", x = pop(regional_obj)) # combine VIU
pop(regional_obj) <- gsub(pattern = "PEN|FRA|JPN", replacement = "JPN", x = pop(regional_obj))              # combine JPN lineage
unique(pop(regional_obj))

pa <- private_alleles(gid = regional_obj)
write.csv(x = pa, file = "03_results/private_alleles.csv", quote = F)

save.image(file = "03_results/post-filters_genind_and_enviro.RData")

# Confirm pops are as expected in obj
as.data.frame(cbind(indNames(obj), as.character(pop(obj))))



# Additional analysis of the per loc stats related to the reason for selecting markers
perloc.FN <- "./03_results/per_locus_stats_2023-08-16.txt"
marker_reason.FN <- "~/Documents/00_sbio/GBMF_UBC_Pacific_oyster/amplicon_panel/ms_oyster_panel_used_for_design_2022-02-07/04_extract_loci/selected_mnames.csv"

perloc.df <- read.table(file = perloc.FN, header = T, sep = "\t")
head(perloc.df)

reason.df <- read.table(file = marker_reason.FN, header = F, sep = ",")
colnames(reason.df) <- c("mname", "reason")
reason.df <- as.data.frame(reason.df)
head(reason.df)

reason.df <- separate(data = reason.df, col = "mname"
                      , into = c("mname", "pos", "nucl"), sep = "_"
                      , remove = TRUE
)
head(reason.df)
reason.df <- reason.df[,c("mname", "reason")]
head(reason.df)

# Combine
all.df <- merge(x = perloc.df, y = reason.df, by = "mname", all.x = T)
head(all.df)

table(all.df[all.df$Hobs > 0.25, "reason"])
length(all.df[all.df$Hobs > 0.25, "reason"])

table(all.df[all.df$Hobs <= 0.25, "reason"])
length(all.df[all.df$Hobs <= 0.25, "reason"])
