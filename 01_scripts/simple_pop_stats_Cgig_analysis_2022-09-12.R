# Main genetic analysis script for the pilot study
# Sutherland Bioinformatics
# Initialized 2022-09-18

# Source simple_pop_stats and choose Pacific oyster

#### 01. Load Data ####
# Loads a genepop to a genind obj
#load_genepop(datatype = "SNP") 
# your data is now obj

# If you are using data consolidated from multiple technical replicates and runs, use '01_scripts/use_multiple_run_data.R' here, then bring the output obj back
#obj.best
#obj <- obj.best

#### 02. Prepare Data ####
## Create a list of individuals for manual addition of population
indiv <- indNames(obj)
indiv.df <- as.data.frame(indiv)
head(indiv.df)

# Separate constituent parts of indiv ID, as these are no longer needed (run, barcode, sample)
indiv.df <- separate(data = indiv.df, col = "indiv", into = c("run", "barcode", "indiv"), sep = "__", remove = T)
head(indiv.df)

indNames(obj) <- indiv.df$indiv

# # Confirm
# test.df <- as.data.frame(indNames(obj.best))
# test.df <- separate(data = test.df, col = "indNames(obj.best)", into = c("run", "barcode", "indiv"), sep = "__", remove = T)
# table(test.df$indiv==indNames(obj))

# How many from each run? 
table(indiv.df$run)

# Keep only necessary column
indiv.df <- indiv.df[, "indiv"]
indiv.df <- as.data.frame(indiv.df)
colnames(indiv.df) <- "indiv"
head(indiv.df)

# Add dummy column to fill manually
indiv.df$pop <- NA

write.table(x = indiv.df, file = "02_input_data/my_data_ind-to-pop.txt"
            , sep = "\t", col.names = T, row.names = F
            , quote = F
            )

# Manually annotate the output file above and populate it with your pop names (no spaces)
# Save as below (tab-delimited) and load
indiv_annot.df <- read.table(file = "02_input_data/my_data_ind-to-pop_annot.txt"
                           , header = T, sep = "\t"
                           #, quote = F
                           )

## Update population names
# Merge with the population annotation, do not sort
indiv_annot_in_order.df <- merge(x = indiv.df, indiv_annot.df, by = "indiv"
                                 , all.x = T, sort = FALSE # very necessary line
                                 )

head(indiv_annot_in_order.df)

# Observe order remained same as (x) above
head(cbind(indiv_annot_in_order.df, indiv.df), n = 10)
tail(cbind(indiv_annot_in_order.df, indiv.df), n = 10)

pop(obj) <- indiv_annot_in_order.df$pop.y
unique(pop(obj))

characterize_genepop(obj)

# ## ASIDE ##
# # Write a little test to be sure
# test <- cbind(indNames(obj), indiv_annot_in_order.df$indiv)
# table(test[,1] == test[ ,2])
# 
# # test <- cbind(indNames(obj), sort(indiv_annot_in_order.df$indiv))
# # table(test[,1] == test[ ,2])
# ## /END/ ASIDE ##

## Population colours
pops_in_genepop <- unique(pop(obj))
pops_in_genepop.df <- as.data.frame(pops_in_genepop)

## Download colours file from previous git repo
# url = "https://raw.githubusercontent.com/bensutherland/ms_oyster_popgen/master/00_archive/my_cols.csv" # only need to run once
destfile <- "00_archive/my_cols.csv"
# download.file(url, destfile)   # only need to run once
my_colours <- read.csv(destfile)
new_pop_colours <- matrix(c("VIU_offspring", "VIU_parent", "red", "pink"), nrow = 2, ncol = 2)
colnames(new_pop_colours) <- c("my.pops", "my.cols")
my_colours <- rbind(my_colours, new_pop_colours)

# Connect colours to empirical populations
colours <- merge(x = pops_in_genepop.df, y =  my_colours, by.x = "pops_in_genepop", by.y = "my.pops"
                 #, sort = F
                 , all.x = T
                 )
colours


#### 03. Characterize missing data (indiv and loci) and filter ####
##### 03.1 Individuals - missing data #####
# # Find missing data per individual
# obj.df <- genind2df(obj)
# #str(obj.df)
# dim(obj.df)
# obj.df[1:5,1:5]
# obj.df[365:370,585:591]
# 
# # Include vector to dataframe of percent missing by individual
# obj.df$ind.per.missing <- NA
# 
# for(i in 1:(nrow(obj.df))){
#   
#   # percent missing             sum of NAs for this row divided by the number of columns, minus two (i.e. pop, per.missing) #TODO: need better way
#   obj.df[i,"ind.per.missing"] <- ( sum(is.na(obj.df[i,])) / (ncol(obj.df) - 2) )
#   #obj.df$ind.per.missing[i] <- ( sum(is.na(obj.df[i,])) / (ncol(obj.df) - 2) )
#   
# }
# 
# dim(obj.df)
# obj.df[1:5, c(1:2, 585:592)]
# obj.df[365:370, c(1:2, 585:592)]

# Combine colours to dataframe for plotting, don't sort
colours
plot_cols.df <- merge(x = obj.df, y = colours, by.x = "pop", by.y = "pops_in_genepop", all.x = T
                      , sort = F
                      )
cols <- plot_cols.df$my.cols

# Plot missing data by individual
pdf(file = "03_results/missing_data_percent_by_ind.pdf", width = 8, height = 5)
plot(obj.df$ind.per.missing, ylab = "Percent missing by individual"
     , col = cols
     , las = 1
     )

abline(h = 0.3, lty = 3)

legend("topright", legend = unique(plot_cols.df$pop)
       , fill = unique(plot_cols.df$my.cols)
       , cex = 0.7
       )
dev.off()


## Filter individuals
# Keep inds with 70% genotyping rate (% missing < 0.3)
keep <- rownames(obj.df[obj.df$ind.per.missing < 0.3, ])

length(keep)
nInd(obj)

obj.filt <- obj[(keep)]
obj.filt

##### 03.2 Loci - missing data #####
# Filter loci based on missing data
obj.df <- genind2df(obj.filt)
obj.df[1:5,1:5]
obj.df <- t(obj.df)
obj.df[1:5,1:5]
obj.df <- obj.df[2:nrow(obj.df),] # remove pop row
obj.df[1:5,1:5]
dim(obj.df)
str(obj.df)

obj.df <- as.data.frame(obj.df)
dim(obj.df)
str(obj.df)
obj.df[1:5,1:5]
obj.df[585:590,290:295]

# Add collector col
obj.df$marker.per.missing <- NA

for(i in 1:(nrow(obj.df))){
  
  # Per marker                      sum all NAs for th emarker, divide by total number markers (#TODO: Confirm or find better method)
  obj.df$marker.per.missing[i] <- ( sum(is.na(obj.df[i,])) / (ncol(obj.df)) )
  
}


# Plot
pdf(file = "03_results/missing_data_by_marker.pdf", width = 5, height = 4)
plot(obj.df$marker.per.missing, xlab = "Marker index", ylab = "Proportion missing data", las = 1)
abline(h = 0.3
       #, col = "grey60"
       , lty = 3, )
dev.off()

# Filter
keep <- rownames(obj.df[obj.df$marker.per.missing < 0.3, ])

# How many loci were removed? 
nLoc(obj.filt) - length(keep)

# Drop loci from genind
obj.all.filt <- obj.filt[, loc=keep]

# Rename back to obj
obj <- obj.all.filt


##### 03.3 Post-missing data filter #####
obj

## View the ind or loc names
inds <- indNames(obj)
loci <- locNames(obj)

# Save out which individuals have passed the filters
write.table(x = inds, file = "03_results/retained_individuals.txt", sep = "\t", quote = F
            , row.names = F, col.names = F
            )

write.table(x = loci, file = "03_results/retained_loci.txt", sep = "\t", quote = F
            , row.names = F, col.names = F
)

##### 03.4 Drop monomorphic loci #####
drop_loci(drop_monomorphic = TRUE)
obj <- obj_filt


##### 03.5 per marker stats and filters #####
## Per locus statistics
per_locus_stats(data = obj)
# note: this writes out to "per_locus_stats*.txt"

# # per-locus heterozygosity (obs)
# obj.summary <- summary(obj)
# str(obj.summary)
# 
# pdf(file = "per_locus_Hobs.pdf", width = 6, height = 5)
# plot(obj.summary[["Hobs"]]
#      , xlab = "Marker (index)"
#      , ylab = "Observed Heterozygosity"
#        )
# abline(h = 0.5, lty = 3)
# dev.off()
# 
# # Which markers are greater than 0.5 heterozygosity? 
# hobs.outliers <- names(obj.summary[["Hobs"]][obj.summary[["Hobs"]] > 0.5])
# hobs.outliers
# 
# keep <- setdiff(x = locNames(obj), y = hobs.outliers)
# 
# # Drop Hobs > 0.5 loci from genind
# obj <- obj[, loc=keep]
# obj

## Hardy-weinberg
hwe_eval(data = obj, alpha = 0.01)
# writes out as HWE_result_alpha_0.01.txt
#TODO: The above Hobs and HWE filters have not been applied to the data yet 


##### 03.5 Post-all filters #####
characterize_genepop(df = obj, N = 30)

# Save out colours
colours
colnames(x = colours) <- c("collection", "colour")
write.csv(x = colours, file = "00_archive/formatted_cols.csv", quote = F, row.names = F)


#### 04. Analysis ####
## Multivariate
# PCA from genind
pca_from_genind(data = obj, PCs_ret = 4, colour_file = "00_archive/formatted_cols.csv")

dapc_from_genind(data = obj, plot_allele_loadings = TRUE, colour_file = "00_archive/formatted_cols.csv")

## Dendrogram
make_tree(bootstrap = TRUE, boot_obj = obj, nboots = 10000, dist_metric = "edwards.dist", separated = FALSE)

## Genetic differentiation
#calculate_FST(format = "genind", dat = obj, separated = FALSE)
calculate_FST(format = "genind", dat = obj, separated = FALSE, bootstrap = TRUE)

## Private alleles
regional_obj <- obj

# Group common pops to get regional private alleles
unique(pop(regional_obj))
pop(regional_obj) <- gsub(pattern = "VIU_offspring|VIU_parent", replacement = "VIU", x = pop(regional_obj))
pop(regional_obj) <- gsub(pattern = "PEN|FRA|JPN", replacement = "JPN", x = pop(regional_obj))
unique(pop(regional_obj))

pa <- private_alleles(gid = regional_obj)
write.csv(x = pa, file = "03_results/private_alleles.csv", quote = F)


#### Convert genepop to Rubias format
# Need to create a stock code file, in the form of
# in the tab-delim format of: 
#collection	repunit
#12Mile_Creek	GoA

stock_code.df <- as.data.frame(unique(pop(obj)))
colnames(stock_code.df) <- "collection"
stock_code.df$repunit <- stock_code.df$collection
stock_code.df
write_delim(x = stock_code.df, file = "00_archive/stock_code.txt", delim = "\t", col_names = T)
micro_stock_code.FN <- "00_archive/stock_code.txt"
# this is for annotate_rubias(), for an unknown reason it requires the name micro_stock_code.FN


sample_to_pop_interp.FN <- "02_input_data/my_data_ind-to-pop_annot.txt"


### TODO: test out the code updates below
genepop_to_rubias_SNP()
# May need some adapting because of pop names and indiv names diff from expected for the pipeline
# Completed 2022-09-19

full_sim(rubias_base.FN = "03_results/rubias_output_SNP.txt", num_sim_indiv = 200, sim_reps = 100)

# Estimating inbreeding (from adegenet tutorial)
obj_PEN <- seppop(x = obj)$PEN
obj_VIU_parent <- seppop(x = obj)$VIU_parent
obj_VIU_offspring <- seppop(x = obj)$VIU_offspring
obj_DPB <- seppop(x = obj)$DPB

# compute the mean inbreeding for each individual and plot
#temp <- inbreeding(x = obj_PEN, N = 100)
#temp <- inbreeding(x = obj_VIU_parent, N = 100)
#temp <- inbreeding(x = obj_VIU_offspring, N = 100)
temp <- inbreeding(x = obj_DPB, N = 100)

class(temp)
head(names(temp))
temp[[1]] # temp is a list of values sampled from the likelihood distribution of each individual; means values are obtained for all indiv using sapply
Fbar <- sapply(temp, mean)
hist(Fbar, col = "firebrick", main = "Average inbreeding in Pendrell")
hist(Fbar, col = "firebrick", main = "Average inbreeding in VIU parents")
hist(Fbar, col = "firebrick", main = "Average inbreeding in VIU offspring")
hist(Fbar, col = "firebrick", main = "Average inbreeding in DPB")


## Per sample heterozygosity

# The following would need extensive coding to make happen
#rubias_to_vcf() # write out, then use instructions here to get per individual heterozygosity in vcftools
# https://github.com/bensutherland/ms_oyster_popgen/blob/master/01_scripts/heterozygosity.sh
# per population heterozygosity





# related would be good to run after here

