### Description: Statistics and Graphing 16S dataset      ###
### Author: Anna Li, date: 17092025                       ###
### Email: anna.li@adelaide.edu.au                        ###

## Load Libraries
library(ape)              # Required for read.tree
library(phyloseq)         # Required for convert relative abundance to phyloseq object
library(vegan)            # Required for adonis2 - beta diversity
library(ggplot2)          # Required for graphing
library(ggpubr)           # Required for stat_compare_means when graphing shannon
library(gridExtra)        # Required for arranging plots (grid.arrange)

## Set up folders
project_dir <- here::here()
results_dir <- file.path(project_dir, "Results")

#====Load phyloseq objects====

physeq_raw <- readRDS(file.path(project_dir, "phyloseq_objects", "physeq_raw.rds"))
physeq_prop <- readRDS(file.path(project_dir, "phyloseq_objects", "physeq_prop.rds"))
physeq_raw_tree <- readRDS(file.path(project_dir, "phyloseq_objects", "physeq_raw_tree.rds"))
metadata_keep <-readRDS(file.path(project_dir, "phyloseq_objects", "metadata_keep.rds"))

#====Estimate alpha and beta diversity in phyloseq and create PC_data.df====
#phyloseq function to check overall richness of data in different measures
plot_richness(physeq_raw)

#alpha diversity in Phyloseq
alpha_div <- estimate_richness(physeq_raw, measures = c("Observed", "Shannon", "Simpson"))
head(alpha_div)
library(dplyr) #for data manipulation/component of tidyverse -  for cleanup of alpha_div
rownames(alpha_div) <- sub("^X", "", rownames(alpha_div))
rownames(alpha_div) <- gsub("\\.", "-", rownames(alpha_div))  # replace dots with dashes

# Beta diversity in PhyloSeq
#bray-curtis
braycurtis_dist = phyloseq::distance(physeq_prop,method="bray") #input physeq_prop for bray curtis analysis as it calculates using composition
braycurtis_pc = pcoa(braycurtis_dist) #make sure to run ape package for pcoa function

#Unifrac
unifrac_dist = phyloseq::distance(physeq_raw_tree,method="unifrac")
unifrac_pc = pcoa(unifrac_dist)

#weighted Unifrac
wunifrac_dist = phyloseq::distance(physeq_raw_tree,method="wunifrac")
wunifrac_pc = pcoa(wunifrac_dist)

#Using Yen Kai's approach to merging data and metadata:
#Creating a single dataframe/table to reference for graphing
PC_data.df = data.frame(metadata_keep,alpha_div,
                        "BC1" = braycurtis_pc$vectors[,1],
                        "BC2" = braycurtis_pc$vectors[,2],
                        "UF1" = unifrac_pc$vectors[,1],
                        "UF2" = unifrac_pc$vectors[,2],
                        "WU1" = wunifrac_pc$vectors[,1],
                        "WU2" = wunifrac_pc$vectors[,2])
group_cols <- c(
  "Donor match FMT"     = "#1b9e77",   # green
  "AutoFMT"             = "#d95f02",   # orange
  "C57/BL6 microbiome"  = "#7570b3",   # purple
  "Saline"              = "#e7298a",   # pink
  "LPJ microbiome"      = "#66a61e"    # olive green
)

# Define the order of the timepoints
#timepoint_levels <- c("S0", "D0", "D15", "F7", "F14", "F21", "Cull", "cecal")

# Add an ordered factor column 
PC_data.df$Timepoint_ordered <- factor(PC_data.df$Timepoint,
                                       #levels = timepoint_levels,
                                       ordered = TRUE)

#Define C57/BL6 and LPJ microbiome as 'timepoint_ordered'
PC_data.df$Timepoint_ordered <- as.character(PC_data.df$Timepoint_ordered)
PC_data.df$Timepoint_ordered[PC_data.df$Sample_ID %in% c("25-03302")] <- "C57/BL6 microbiome"
PC_data.df$Timepoint_ordered[PC_data.df$Sample_ID %in% c("25-03303")] <- "LPJ microbiome"
PC_data.df$Timepoint_ordered <- factor(PC_data.df$Timepoint_ordered, 
                                       levels = c("C57/BL6 microbiome","S0","D0","D15","F7","F14","F21","Cull","cecal","LPJ microbiome"))
#Can omit cecal here if don't  want to include in longitudinal graphs (Timepoint_ordered arguments)

#====PERMANOVA tests=====
#Formula syntax: adonis2(distance_matrix ~ factor_variable, data = metadata)
#For more permutations: adonis2(braycurtis_dist ~ Group, data = metadata, permutations = 9999)
#To test multiple variables: adonis2(braycurtis_dist ~ Group + Timepoint, data = metadata)

#braycurtis
adonis2(braycurtis_dist ~ Timepoint,metadata_keep) #sig = ***
adonis2(braycurtis_dist ~ Group,metadata_keep) #sig = ***
adonis2(braycurtis_dist ~ Timepoint + Group, data = metadata_keep) #sig = ***

#Additional analysis I can do: pairwiseAdonis2 package to look at which groups differ

#unifrac
adonis2(unifrac_dist ~ Timepoint,metadata_keep)

#weighted unifrac
adonis2(wunifrac_dist ~ Timepoint,metadata_keep)

##=====Graphing alpha diversity====
#Create object without FMT samples for alpha diversity plots
Recipient_data.df <- PC_data.df[!rownames(PC_data.df)%in% c("25-03302", "25-03303"), ]

#----Graphing shannon diversity across timepoints----
Shannon_timeline_plot <- ggplot(Recipient_data.df, aes(x = Timepoint_ordered, y = Shannon, fill = Timepoint_ordered)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_line(aes(group = Mouse_ID), alpha = 0.2, linewidth = 0.2) +
  geom_point(aes(shape = Death_timepoint, colour = Group), size = 2.5, alpha = 0.9) +
  theme_minimal(base_size = 14) +
  labs(title = "Shannon Diversity per Timepoint",
       x = "Timepoint",
       y = "Shannon Diversity Index") +
  theme(legend.position = "right")

Shannon_timeline_plot

#----Plotting Shannon's diversity - by group at specific timepoints----
# Generate all pairwise combinations of groups/define pairwise comparisons
comparisons <- combn(unique(Recipient_data.df$Group), 2, simplify = FALSE)

#Plotting Shannon's diversity - at D15 by group
# Subset only D15 samples
alpha_div_D15 <- Recipient_data.df %>%
  filter(Timepoint_ordered == "D15")
# Kruskal-Wallis test across all groups at cull
kruskal_result <- kruskal.test(Shannon ~ Group, data = alpha_div_D15)
print(kruskal_result)

# Boxplot with kruskal-Wallis comparisons and wilcox (can run separately). 
shannon_D15_plot <- ggplot(alpha_div_D15, aes(x = Group, y = Shannon, fill = Group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA,) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7) +
  theme_minimal(base_size = 14) +
  labs(title = "Shannon Diversity at D15 by Group",
       x = "",
       y = "Shannon Diversity Index") +
  theme(legend.position = "none") +
  stat_compare_means(method = "kruskal.test", label.y = max(alpha_div_D15$Shannon) + 0.2) + # global test
  stat_compare_means(method = "wilcox.test", 
                     comparisons = combn(unique(alpha_div_D15$Group), 2, simplify = FALSE),
                     label = "p.signif")
shannon_D15_plot

#----plotting Shannon's diversity - by group for all individual timepoints---- 

shannon_all_plot <- ggplot(Recipient_data.df, aes(x = Group, y = Shannon, fill = Group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(#aes(color = HSCT_challenge), #to colour the points by metadata you like
              width = 0.2, size = 2, alpha = 0.7) +   # color by metadata input i.e. Study_ID
  facet_wrap(~ Timepoint_ordered) +   # cycle through timepoints
  theme_minimal(base_size = 14) +
  labs(title = "Shannon Diversity Across Timepoints by Treatment Group",
       x = "",
       y = "Shannon Diversity Index") +
  #scale_color_manual(values = c("Mild" = "yellow", "Moderate" = "orange", "Robust" = "purple")) +   #choose the metadata names and colour
  theme(legend.position = "none") +
  stat_compare_means(method = "kruskal.test", 
                     label.y.npc = "top", label.x.npc = "left") + 
  stat_compare_means(method = "wilcox.test",
                     comparisons = combn(unique(Recipient_data.df$Group), 2, simplify = FALSE),
                     label = "p.signif")

shannon_all_plot

#----plotting Shannon's diversity as a raincloud----
library(tidyverse)
library(ggdist)
library(gghalves)

#stats
#Creating mean statistics for groups across timepoints 
library(dplyr)

shannon_summary <- Recipient_data.df %>%
  group_by(Group, Timepoint_ordered) %>%
  summarise(
    mean_shannon = mean(Shannon, na.rm = TRUE),
    se_shannon   = sd(Shannon, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

#plot
raincloud_shannon <- ggplot(Recipient_data.df, 
                            aes(x = Timepoint_ordered, y = Shannon, fill = Timepoint_ordered)) +
  # Half-violin (density)
  stat_halfeye(
    adjust = 0.7, width = 0.7, .width = 0, justification = -0.2,
    point_interval = median_qi,
    alpha = 0.6
  ) +
  # Jittered points
  geom_point(
    aes(color = Group), 
    position = position_jitter(width = 0.1, height = 0),
    alpha = 0.6, size = 2
  ) +
  # Boxplot (optional; remove if you want summary dot only)
  geom_boxplot(
   width = 0.15, outlier.shape = NA,
   alpha = 0.7, color = "black"
  ) +
  # Mean Shannon line per group
  geom_line(data = shannon_summary,
            aes(x = Timepoint_ordered, y = mean_shannon, group = Group, color = Group),
            linewidth = 1.2) +
  
  geom_point(data = shannon_summary,
             aes(x = Timepoint_ordered, y = mean_shannon, color = Group),
             size = 3) +
  theme_minimal(base_size = 14) +
  labs(title = "Shannon Diversity by Timepoint – Raincloud Plot",
       x = "",
       y = "Shannon Diversity") +
  theme(legend.position = "right")

#plot
raincloud_shannon

#====Graphing beta diversity====
#----BC_group_plot - all samples coloured by group----
BC_group_plot <- ggplot(PC_data.df, aes(x = BC1, y = BC2, colour = Group, size = Group)) +
  geom_point(alpha = 0.8) +
  scale_size_manual(values = c(
    "AutoFMT"             = 2,
    "Donor match FMT"     = 2,
    "Saline"              = 2,
    "C57/BL6 microbiome"  = 6,
    "LPJ microbiome"      = 6
  )) +
  scale_color_manual(values = group_cols) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Bray–Curtis PCoA",
    x = "PCoA Axis 1",
    y = "PCoA Axis 2"
  )

BC_group_plot

#----BC_timepoint_plot - all samples coloured by timepoint (not group)----

BC_timepoint_plot <- ggplot(PC_data.df, aes(x = BC1, y = BC2, 
                                  colour = Timepoint_ordered, 
                                  size = Group)) + #Inside aes() you don’t set exact colors/sizes; instead, you tell ggplot which variable controls them
  geom_point(alpha = 0.8) +
  scale_size_manual(values = c(
    "AutoFMT"             = 2,
    "Donor match FMT"     = 2,
    "Saline"              = 2,
    "C57/BL6 microbiome"  = 6,
    "LPJ microbiome"      = 6
  ), guide = "none") +  # hide size legend if not needed
  scale_color_viridis_d(option = "plasma") +  # discrete palette for ordered factors such as (Timepoint_ordered)
  theme_minimal(base_size = 14) +
  labs(
    title = "Bray–Curtis PCoA (coloured by Timepoint)",
    x = "PCoA Axis 1",
    y = "PCoA Axis 2",
    colour = "Timepoint"
  )

BC_timepoint_plot

#----BC_PLOTS - graphs for each timepoint----

#chatgpt helped here...
special_samples <- c("25-03302", "25-03303")
# Keep only real timepoints (for plotting individual graphs)
timepoints <- setdiff(unique(PC_data.df$Timepoint_ordered),
                      c("LPJ microbiome", "C57/BL6 microbiome"))

timepoints <- setdiff(levels(PC_data.df$Timepoint_ordered),
                      c("LPJ microbiome", "C57/BL6 microbiome"))

BC_PLOTS <- list()

for(tp in timepoints){
  df_tp <- PC_data.df %>% filter(Timepoint_ordered == tp)
  
  BC_PLOTS[[tp]] <- ggplot(df_tp, aes(x = BC1, y = BC2, colour = Group, size = Group)) +
    geom_point(alpha = 0.8, stroke = 0.5) +
    # ellipses for each group at that timepoint
    stat_ellipse(level = 0.95, size = 0.5) +
    # overlay special samples: keep group color, but larger + black outline
    geom_point(data = PC_data.df %>% filter(Sample_ID %in% special_samples),
               aes(x = BC1, y = BC2, colour = Group), 
               size = 5, shape = 21, stroke = 1.5, fill = NA) +
    scale_color_manual(values = group_cols) +
    scale_size_manual(values = c(
      "AutoFMT"            = 2,
      "Donor match FMT"    = 2,
      "Saline"             = 2,
      "C57/BL6 microbiome" = 6,
      "LPJ microbiome"     = 6
    ), guide = "none") +
    theme_minimal(base_size = 14) +
    labs(title = paste("Bray–Curtis PCoA:", tp),
         x = "PCoA Axis 1",
         y = "PCoA Axis 2")
}
BC_PLOTS

library(cowplot)
BC_time_plots <- plot_grid(plotlist = BC_PLOTS[timepoints],
                           ncol = 2, labels = "AUTO")
BC_time_plots

#====save graphs for export====
#Create consistent filepath and parameters for saving graphs

save_plot <- function(plot_obj, filename, width = 8, height = 6) {
  ggsave(
    filename = file.path(results_dir, filename),
    plot = plot_obj,
    width = width,
    height = height
  )
}

save_plot(Shannon_timeline_plot, "Shannon_timeline_plot.pdf")
save_plot(BC_time_plots,"BrayCurtis_PCoA_allTimepoints.pdf", width = 12, height = 10)

#====save the current workspace====
save.image(file = "GvHDfmt5_stats-plots.RData")
saveRDS(PC_data.df, file=file.path(project_dir, "Phyloseq_objects","PC_data.df.rds"))

