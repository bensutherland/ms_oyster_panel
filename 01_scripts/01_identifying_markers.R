# Script to identify top markers for Pacific oyster panel design
# B. Sutherland (Sutherland Bioinformatics)
# 2022-01-07

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
input.FN <- "02_input_data//populations_single_snp_HWE.raw"
nmarkers_to_keep <- 300


#### 1. Import data ####
print(paste0("Loading data from ", input.FN))
my.data <- read.PLINK(file = input.FN)
my.data


#### 2. Observe data (allele frequencies) ####
# Plot instances of minor allele across individuals and loci
png(file = paste0(output.dir, "maf_glPlot_all.png"), width = 924, height = 600)
glPlot(x = my.data, posi="topleft") 
dev.off()
# note: in this data, it appears that some markers have the second allele as the more common allele

# Create density plot of minor allele frequencies
pdf(file = paste0(output.dir, "maf_hist_all.pdf"), width = 6, height = 4)
myFreq <- glMean(my.data)
hist(myFreq
     , proba=T # means not frequency based, but rather probability densities
     , col="gold", xlab = "Allele frequencies"
     , main = ""
     , ylim = c(0,50)
     , ylab = "Density of second allele frequencies"
)
text(x = 0.4, y = 7, labels = paste(nLoc(my.data), " loci", sep = "" ))
temp <- density(myFreq)
lines(temp$x, temp$y, lwd=3)
dev.off()


####3. Convert genlight into the essential formats ####
# Convert genlight to matrix
my.data.mat <- as.matrix(my.data)
my.data.mat[1:5,1:5]

# Translate the number of minor allele to genind format
my.data.mat[my.data.mat == 0] <- "1/1" #homozygote reference
my.data.mat[my.data.mat == 1] <- "1/2" #heterozygote
my.data.mat[my.data.mat == 2] <- "2/2" #homozygote alternate
my.data.mat[1:5,1:5]

# Convert matrix to genind
my.data.gid <- df2genind(my.data.mat, sep = "/", ploidy = 2) # convert df to genind

# Transfer pop attributes
pop(my.data.gid) <- pop(my.data) 

# Show sample size per population
table(pop(my.data.gid))
pdf(file = paste0(output.dir, "sample_size_barplot.pdf"), width = 8, height = 5)
par(mfrow=c(1,1), mar=c(8,5,3,3))
barplot(table(pop(my.data.gid)), col=funky(17)
        #, las=3, las = 1
        , las=2
        , xlab=""
        , ylab="Sample size"
        , ylim = c(0,40))
#abline(h = c(10,20,30), lty=2)
dev.off()

# Convert genind to hierfstat
all.data.hf <- genind2hierfstat(my.data.gid)
rownames(all.data.hf) <- indNames(my.data.gid)

# Convert genind to genepop
all_data.gen <- genind2genpop(x = my.data.gid, pop = pop(my.data.gid))

## Summary Now have data in either hf format, genind, and genepop
## all.data.hf # hf
## my.data.gid # genind
## all_data.gen # genepop


#### 4. Filter to retain specific populations (BC naturalized) ####
sep.obj <- seppop(x = my.data.gid)
names(sep.obj)

# Create list of various combinations of repooled pops
target_pops.gid <- repool(sep.obj$HIS, sep.obj$PEN, sep.obj$PIP, sep.obj$SER)
# rm(sep.obj)


#### 5. MAF statistics and filter #####
# Calculate MAF
maf.vec <- minorAllele(target_pops.gid)
length(maf.vec) # how many markers present? 

# Construct MAF into dataframe
maf.df <- as.data.frame(maf.vec)
maf.df$mname <- rownames(maf.df)
maf.df <- maf.df[,c("mname","maf.vec")]
head(maf.df, n = 50)

# Note: MAF in some cases includes markers that have 1.0 MAF, which should also be removed (alt being called main allele?)

# Filter based on MAF
# Which markers have global MAF < 0.01? 
markers_to_remove.df <- maf.df[maf.df$maf.vec < 0.01, ]
head(markers_to_remove.df)
nrow(markers_to_remove.df)

# Which markers have global MAF > 0.99? 
markers_to_remove_alt.df <- maf.df[maf.df$maf.vec > 0.99, ]
head(markers_to_remove_alt.df)
nrow(markers_to_remove_alt.df)

# Create dataframe with markers to drop
markers_to_remove.df <- rbind(markers_to_remove.df, markers_to_remove_alt.df)

# Note: could technically introduce drop_monomorphic here

# Reporting the number markers being dropped
print(
  paste0("If correct, will remove ", nrow(markers_to_remove.df) 
         ,  " markers from the input of "
         , nLoc(target_pops.gid)
         , " and will retain a total of " 
         , nLoc(target_pops.gid) - nrow(markers_to_remove.df)
         , " markers"
  )
  )
  
# Find the names of the markers to retain
mnames_to_retain.vec <- setdiff(x = locNames(target_pops.gid), y = markers_to_remove.df$mname)
length(mnames_to_retain.vec)

# Use retain list to keep the MAF-filtered markers
target_pops_maf_filt.gid <- target_pops.gid[loc=mnames_to_retain.vec]
target_pops_maf_filt.gid

#### 6. Heterozygosity statistics ####
# Calculate descriptive statistics on genind file
div <- summary(target_pops_maf_filt.gid) # genetic diversity w/ adegenet
str(div) # What is contained within this summary object
# contains: sample size, sample size by pop, number loci and names, population IDs, % NA, Hobs, Hexp

# Create dataframe of observed heterozygosity
hobs.df <- data.frame(div[["Hobs"]])
dim(hobs.df)
colnames(hobs.df) <- "Hobs"
hobs.df$mname <- rownames(hobs.df)
hobs.df <- hobs.df[,c(2,1)]
hobs.df <- hobs.df[order(hobs.df$Hobs, decreasing = TRUE), ] #sort
head(hobs.df)

# Create frequency plot of Hobs
pdf(file = paste0(output.dir, "hobs_hist.pdf"), width = 10, height = 5)
Hobs.freq <- div[["Hobs"]]
hist(Hobs.freq
     #, proba=T # means not frequency based, but rather probability densities
     , freq = TRUE
     , col="grey"
     , xlab = "Hobs"
     , main = ""
     #, ylim = c(0,0.8)
     , ylab = "Freq. of Hobs (# markers)"
     , breaks = 50
     , las = 1
)

# Note: assumes sorted, add a line indicating the position of the top X markers 
# The lower limit of the top n markers under 0.5 hobs
abline(v = tail(head(hobs.df[hobs.df$Hobs < 0.5,], n = nmarkers_to_keep), n = 1)[2], lty = 2)
text(x = 0.6, y = 1500
     , labels = paste0("Top 300 markers "
                       , '\n'
                       , round(x = tail(head(hobs.df[hobs.df$Hobs < 0.5,], n = nmarkers_to_keep), n = 1)[2], digits = 2)
                       , " < Hobs < 0.5 "
                       )
)
abline(v = 0.5, lty = 2)
dev.off()

# Add the marker names to the hobs, filter, and only keep the target number
nrow(hobs.df)
hobs_selected.df <- hobs.df[hobs.df$Hobs < 0.5, ] # Remove those above 0.5
nrow(hobs_selected.df)
hobs_selected.df <- head(hobs_selected.df, n = nmarkers_to_keep)
nrow(hobs_selected.df)

# Save out hobs selected markers
# write.csv(x = hobs_selected.df, file = "03_marker_selection/top_hobs_markers.csv", quote = F, row.names = F)
## note: we use the all markers file to identify selected markers, not this file. 


#### 7. Differentiation statistics (Fst) ####
# Calculate Fit, Fst, and Fis per locus
per_locus_stats <- Fst(as.loci(target_pops_maf_filt.gid))
per_locus_stats.bck <- per_locus_stats # create backup

# Format into dataframe in order of descending Fst
per_locus_stats.df <- as.data.frame(per_locus_stats)
str(per_locus_stats)
per_locus_stats.df$mname <- rownames(per_locus_stats.df)
per_locus_stats.df <- per_locus_stats.df[, c("mname", "Fit", "Fst", "Fis")]
head(per_locus_stats.df)
per_locus_stats.df <- per_locus_stats.df[order(per_locus_stats.df$Fst, decreasing = TRUE), ] #sort
head(per_locus_stats.df)
dim(per_locus_stats.df)

# Write out all per locus stats
#write.csv(x = per_locus_stats.df, file = "03_marker_selection/all_markers_per_locus_stats.csv", quote = F, row.names = F)
## note: we use the merged file below, not this file, for marker selection

# Create density plot of Fst
pdf(file = paste0(output.dir, "fst_hist.pdf"), width = 10, height = 5)
fst.freq <- per_locus_stats.df$Fst
hist(fst.freq
     #, proba=T # means not frequency based, but rather probability densities
     , freq = TRUE
     , col="grey"
     , xlab = "Fst"
     , main = ""
     #, ylim = c(0,0.8)
     , ylab = "Freq. of Fst vals (# markers)"
     , breaks = 50
     , las = 1
)

# Assumes sorted, add a line indicating the position of the top X markers 
# The lower limit of the top n markers under 0.5 hobs
abline(v = tail(head(per_locus_stats.df, n = nmarkers_to_keep), n = 1)["Fst"], lty = 2)
text(x = 0.3, y = 3000
     , labels = paste0("Top 300 markers "
                       , '\n'
                       , "Fst > "
                       , round(x = tail(head(per_locus_stats.df, n = nmarkers_to_keep), n = 1)["Fst"], digits = 3)
                       
     )
)

dev.off()

# Stats
# Lowest FST included: 
tail(head(per_locus_stats.df, n = nmarkers_to_keep), n = 1)["Fst"]

# Lowest FST included: 
head(per_locus_stats.df, n = 1)["Fst"]


#### 8. All per-locus statistics together ####
# hobs.df, per_locus_stats.df, maf.df
dim(hobs.df) # maf filtered loci
dim(per_locus_stats.df) # maf filtered loci
dim(maf.df) # all loci

# Merge
all_data.df <- merge(x = per_locus_stats.df, y = maf.df, by = "mname")
dim(all_data.df)
head(all_data.df)

# Merge
all_data.df <- merge(x = all_data.df, y = hobs.df, by = "mname")
dim(all_data.df)
head(all_data.df)

# Write full output
write.csv(x = all_data.df, file = "03_marker_selection/all_markers_per_locus_stats_incl_hobs_MAF.csv", quote = F, row.names = F)
# This is the file that will be used for choosing markers out of the genome 


#### 9. Compare statistics against each other ####
par(mfrow = c(2,2))
#pdf(file = paste0(output.dir, "comparative_stats.pdf"), width = 8, height = 8)

plot(x = all_data.df$Fst, y = all_data.df$maf.vec, xlab = "Fst", ylab = "Minor Allele Frequency")
abline(v = tail(head(per_locus_stats.df, n = nmarkers_to_keep), n = 1)["Fst"], lty = 2)
#dev.off()

#pdf(file = paste0(output.dir, "hobs_vs_fst.pdf"), width = 10, height = 5)
plot(x = all_data.df$Fst, y = all_data.df$Hobs, ylab = "Hobs", xlab = "Fst")
abline(v = tail(head(per_locus_stats.df, n = nmarkers_to_keep), n = 1)["Fst"], lty = 2)
#dev.off()

#pdf(file = paste0(output.dir, "maf_vs_hobs.pdf"), width = 10, height = 5)
plot(x = all_data.df$Hobs, y = all_data.df$maf.vec, xlab = "Observed heterozygosity", ylab = "Minor Allele Frequency")
abline(v = tail(head(hobs.df[hobs.df$Hobs < 0.5,], n = nmarkers_to_keep), n = 1)[2], lty = 2)
abline(v = 0.5, lty = 2)

#dev.off()

#save out manually as 8x8 (can't seem to output as 2x2 automatically, will try again #todo)


#### 10. Write out Rdata
##### 9. Save Results ####
save.image(file = paste0(output.dir, "adegenet_output.RData"))


#### Next Steps: #####
# Now you have both the top FST markers, and the best heterozygosity markers
# The data was already filtered for HWE and MAF (within selected pops)
# The next step is to take the lists of top Fst markers and top Hobs markers, then use these to extract from the vcf and the reference genome

# An additional next step is to go to the Rscript private_alleles.r, which will require you to load the image
#  and use the adegenet_output.RData's my.data.gid to work on some private allele work

