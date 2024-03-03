## Compare HOBS between VIU parentage dataset and the CHR8 parentage dataset
# 2024-02-28
# B. Sutherland

# Load a simple_pop_stats for functions

#### 01. Plot VIU vs. OCP HOBS ####
# Obtain parentage genind from pilot, and calculate per locus HOBS
data_pilot.FN <- "~/Documents/00_sbio/GBMF_UBC_Pacific_oyster/amplicon_panel/pilot_study/simple_pop_stats_v.0.5_2024-02-20/03_results/post-filters_prepared_for_parentage_rubias_built_problem_genos_rem_rubias_built.RData"
load(data_pilot.FN)
parentage_obj_pilot <- obj_filt
per_locus_stats(data = parentage_obj_pilot)
per_loc_stats_pilot.df <- per_loc_stats.df
nrow(per_loc_stats_pilot.df) # 347

# Obtain parentage genind from CHR8, and calculate per locus HOBS
data_OCP23.FN <- "~/Documents/00_sbio/VIU/VIU_oyster/CHR8/OCP23/OCP23_parentage_v.0.3/simple_pop_stats_OCP23_v.0.3/03_results/post-filters_prepared_for_parentage_rubias_built_problem_genos_rem_rubias_built.RData"
load(data_OCP23.FN)
parentage_obj_OCP23 <- obj_filt
per_locus_stats(data = parentage_obj_OCP23) # NOTE: this will overwrite the original per locus stats output in the folder
per_loc_stats_OCP23.df <- per_loc_stats.df
nrow(per_loc_stats_OCP23.df) # 300

# View data
head(per_loc_stats_OCP23.df)

# Combine both datasets
per_loc_stats_both_datasets.df <- merge(x = per_loc_stats_pilot.df, y = per_loc_stats_OCP23.df, by = "mname")
head(per_loc_stats_both_datasets.df)
dim(per_loc_stats_both_datasets.df) # 248 shared
colnames(per_loc_stats_both_datasets.df) <- gsub(pattern = "\\.x", replacement = ".pilot", x = colnames(per_loc_stats_both_datasets.df))
colnames(per_loc_stats_both_datasets.df) <- gsub(pattern = ".y", replacement = ".OCP23", x = colnames(per_loc_stats_both_datasets.df))
head(per_loc_stats_both_datasets.df)


## Plotting
#pdf(file = "03_results/parentage_dataset_comparative_hobs.pdf", width = 5, height = 5)
p <- ggplot(data = per_loc_stats_both_datasets.df, aes(x = Hobs.pilot, y = Hobs.OCP23)) +
  geom_point() +
  theme_bw() +
  xlab(bquote("per locus "*H[OBS]*" (VIU OFR5+OFR1)")) + 
  ylab(bquote("per locus "*H[OBS]*" (MBP CHR8)")) +
  xlim(0, 0.9) + 
  ylim(0, 0.9)
p

num_loci <- paste0(nrow(per_loc_stats_both_datasets.df), " loci")

p <- p + geom_text(x = 0.1, y = 0.77, label = num_loci, check_overlap = T)
p

# rename to keep
hobs_pilot_v_OCP23.plot <- p
hobs_pilot_v_OCP23.plot


#### 02. Plot pilot HOBS vs. FST ####
data_pilot_popgen.FN <- "03_results/post-filters_prepared_for_parentage.RData"
load(data_pilot_popgen.FN)
# Get back original, all data
obj <- repool(separated_pops$VIU_F2, separated_pops$VIU_F1, separated_pops$VIU_F0, separated_pops$DPB
       , separated_pops$PEN, separated_pops$GUR, separated_pops$FRA, separated_pops$JPN, separated_pops$CHN
          )
obj # should be 425 loci and 312 inds after filters
per_locus_stats(obj)
per_loc_stats_popgen.df <- per_loc_stats.df
nrow(per_loc_stats_popgen.df) # 425 loci
head(per_loc_stats_popgen.df)

## Plotting
p <- ggplot(data = per_loc_stats_popgen.df, aes(x = Fst, y = Hobs)) +
  geom_point() +
  theme_bw() +
  xlab(bquote("per locus "*F[ST]*" (pilot data)")) + 
  ylab(bquote("per locus "*H[OBS]*" (pilot data)")) +
  xlim(0, 0.4) + 
  ylim(0, 0.7)
p

num_loci <- paste0(nrow(per_loc_stats_popgen.df), " loci")

p <- p + geom_text(x = 0.25, y = 0.65, label = num_loci, check_overlap = T)
p

# rename to keep
pilot_hobs_v_fst.plot <- p
pilot_hobs_v_fst.plot


#install.packages(ggpubr)
library(ggpubr)
final.figure <- ggarrange(pilot_hobs_v_fst.plot, hobs_pilot_v_OCP23.plot
                          , labels = c("A", "B"))
final.figure

pdf(file = "03_results/multipanel_per_locus_stats_comparisons.pdf", width = 8.5, height = 4)
print(final.figure)
dev.off()
