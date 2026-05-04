### Description: Statistics and Graphing 16S dataset      ###
### Author: Anna Li, date: 17092025 - 27042026            ###
### Email: anna.li@adelaide.edu.au                        ###

## Load Libraries
library(ape)              # Required for read.tree
library(phyloseq)         # Required for convert relative abundance to phyloseq object
library(vegan)            # Required for adonis2 - beta diversity
library(ggplot2)          # Required for graphing
library(ggpubr)           # Required for stat_compare_means when graphing shannon
library(gridExtra)        # Required for arranging plots (grid.arrange)
library(dplyr)
library(pairwiseAdonis)

## Set up folders
project_dir <- here::here()
results_dir <- file.path(project_dir, "Results")

#----Calculate alpha diversity----
alpha_div <- estimate_richness(physeq_raw, measures = c("Observed", "Shannon", "Simpson"))
head(alpha_div)
library(dplyr) #for data manipulation/component of tidyverse -  for cleanup of alpha_div
rownames(alpha_div) <- sub("^X", "", rownames(alpha_div))
rownames(alpha_div) <- gsub("\\.", "-", rownames(alpha_div))  # replace dots with dashes

#----Calculate Beta diversity----
#bray-curtis
braycurtis_dist = phyloseq::distance(physeq_prop,method="bray") #input physeq_prop for bray curtis analysis as it calculates using composition
braycurtis_pc = pcoa(braycurtis_dist)

#Unifrac
unifrac_dist = phyloseq::distance(physeq_raw_tree,method="unifrac")
unifrac_pc = pcoa(unifrac_dist)

#weighted Unifrac
wunifrac_dist = phyloseq::distance(physeq_raw_tree,method="wunifrac")
wunifrac_pc = pcoa(wunifrac_dist)

#----Create a single dataframe/table to reference for graphing----
PC_data.df = data.frame(map_data,alpha_div,
                        "BC1" = braycurtis_pc$vectors[,1],
                        "BC2" = braycurtis_pc$vectors[,2],
                        "UF1" = unifrac_pc$vectors[,1],
                        "UF2" = unifrac_pc$vectors[,2],
                        "WU1" = wunifrac_pc$vectors[,1],
                        "WU2" = wunifrac_pc$vectors[,2])
group_cols <- c(
  "AutoFMT"             = "#d95f02",   # orange
  "Donor match FMT"     = "#1b9e77",   # green
  "Saline"              = "#0080C0",   # blue
  "C57/BL6 microbiome"  = "#e7298a",   # pink
  "LPJ microbiome"      = "#66a61e"    # olive green
)

## Adding a "timepoint_ordered" column for longitudinal graphs
# Define the order of the timepoints
timepoint_levels <- c("S0", "D0", "D15", "F7", "F14", "F21", "Cull", "Endpoint cecal")

# Add an ordered factor column 
PC_data.df$Timepoint_ordered <- factor(PC_data.df$Timepoint,
                                       levels = timepoint_levels,
                                       ordered = TRUE)

#Define C57/BL6 and LPJ microbiome in 'timepoint_ordered'
PC_data.df$Timepoint_ordered <- as.character(PC_data.df$Timepoint_ordered)
PC_data.df$Timepoint_ordered[PC_data.df$Sample_ID %in% c("25-03302")] <- "C57/BL6 microbiome"
PC_data.df$Timepoint_ordered[PC_data.df$Sample_ID %in% c("25-03303")] <- "LPJ microbiome"
PC_data.df$Timepoint_ordered <- factor(PC_data.df$Timepoint_ordered, 
                                       levels = c("C57/BL6 microbiome","S0","D0","D15","F7","F14","F21","Cull","Endpoint cecal","LPJ microbiome"))
#Can omit cecal here if don't  want to include in longitudinal graphs (Timepoint_ordered arguments)

#=====Graphing alpha diversity====
# Create object without FMT samples for alpha diversity plots
Recipient_data.df <- PC_data.df[!rownames(PC_data.df)%in% c("25-03302", "25-03303"), ]
# Remove "cecal" samples (only plotting stool samples)
Recipient_data.df <- Recipient_data.df %>%
  filter(Timepoint_ordered != "Endpoint cecal")

# Create object where data is only included if Death_timepoint column = endpoint (mouse survived through to end of experiment) and remove cecal
Recipient_data_endpoint.df <- Recipient_data.df %>%
  filter(Death_timepoint == "endpoint", Timepoint_ordered != "Endpoint cecal")

#Calculating mean for groups across timepoints 
shannon_summary <- Recipient_data.df %>%
  group_by(Group, Timepoint_ordered) %>%
  summarise(
    mean_shannon = mean(Shannon, na.rm = TRUE),
    se_shannon   = sd(Shannon, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

#----Graphing shannon diversity across timepoints----
Shannon_timeline_plot <- ggplot(Recipient_data.df, aes(x = Timepoint_ordered, y = Shannon, fill = Timepoint_ordered)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_line(aes(group = Mouse_ID), alpha = 0.2, linewidth = 0.2) +
  geom_point(aes(shape = Death_timepoint, colour = Group), size = 2.5, alpha = 0.9) +
  scale_shape_manual(values = c(
    "endpoint" = 16,            # circle
    "non-endpoint" = 17,        # triangle
    "died preF7" = 15           # square
  )) +
  theme_minimal(base_size = 14) +
  labs(title = "Shannon's diversity per timepoint",
       x = "Timepoint",
       y = "Shannon Diversity Index") +
  labs(subtitle = "All mice throughout experiment")
  theme(legend.position = "right")

Shannon_timeline_plot

#should refer to Feargal's code to see how to add stats across the timepoints

#----Graphing Shannon's diversity - by group at specific timepoints----
# Generate all pairwise combinations of groups/define pairwise comparisons
comparisons <- combn(unique(Recipient_data.df$Group), 2, simplify = FALSE)

#Plotting Shannon's diversity - at D15 by group
# Subset only D15 samples
alpha_div_D15 <- Recipient_data.df %>%
  filter(Timepoint_ordered == "D15")
# Kruskal-Wallis test across all groups at cull
kruskal_result <- kruskal.test(Shannon ~ Group, data = alpha_div_D15)
print(kruskal_result)

# Boxplot with kruskal-Wallis comparisons and wilcox
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

#----Graphing Shannon's diversity - by group for all individual timepoints---- 

Shannon_all_plot <- ggplot(Recipient_data.df, aes(x = Group, y = Shannon, fill = Group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(#aes(color = HSCT_challenge), #to colour the points by any metadata you like
              width = 0.2, size = 2, alpha = 0.7) +   # color by metadata input i.e. Study_ID
  facet_wrap(~ Timepoint_ordered) +   # cycle through timepoints
  theme_minimal(base_size = 14) +
  labs(title = "Shannon's diversity across timepoints by treatment group",
       x = "",
       y = "Shannon Diversity Index") +
  #scale_color_manual(values = c("Mild" = "yellow", "Moderate" = "orange", "Robust" = "purple")) +   #choose the metadata names and colour
  theme(legend.position = "none") +
  stat_compare_means(method = "kruskal.test", 
                     label.y.npc = "top", label.x.npc = "left") + 
  stat_compare_means(method = "wilcox.test",
                     comparisons = combn(unique(Recipient_data.df$Group), 2, simplify = FALSE),
                     label = "p.signif")

Shannon_all_plot
#demonstrates that there are no statistic differences in alpha diversity at each timepoint

#----Graphing Shannon's diversity as a raincloud----
library(tidyverse)
library(ggdist)
library(gghalves)

#Plot: raincloud_shannon
raincloud_shannon_plot <- ggplot(Recipient_data.df, 
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
  labs(subtitle = "All mice throughout experiment") +
  theme(legend.position = "right")

raincloud_shannon_plot

#====Graphing beta diversity====
#----Bray Curtis global statistics / PERMANOVA----
#Note: use recipient_data for stats (don't factor in donor samples for global stats on trajectory yet)
#Align Recipient_data.df with braycurtis_distance 
braycurtis_matrix <- as.matrix(braycurtis_dist)
samples <- intersect(rownames(braycurtis_matrix), rownames(Recipient_data.df))

Recipient_data.df <- Recipient_data.df[samples, ]
braycurtis_subset <- as.dist(braycurtis_matrix[samples, samples])

# Make sure that timepoint_ordered is a factor not a numeric ^should have been coded before but redo here
Recipient_data.df$Timepoint_ordered <- factor(
  Recipient_data.df$Timepoint_ordered,
  levels = c("S0","D0","D15","F7","F14","F21","Cull"),
  ordered = TRUE
)

###Stats: (Question: is variation in gut microbiome explained by time?/Do microbiome compositions differ between timepoints?)
adonis2(braycurtis_subset ~ Timepoint_ordered, 
        data = Recipient_data.df,
        permutations = 999,
        strata = Recipient_data.df$Mouse_ID) #permutations are constrained within each mouse/repeated measures

###Stats: (Question: does autoFMT vs Donor-matchFMT vs Saline lead to different longitudinal changes?)
adonis2(braycurtis_subset ~ Group * Timepoint_ordered,
        data = Recipient_data.df,
        permutations = 999,
        strata = Recipient_data.df$Mouse_ID,
    by = "margin") #term level output

#checking inter dispersion within the time and group (PERMANOVA detects group centroids i.e. biological or dispersion i.e. variances)
betadisper(braycurtis_subset, Recipient_data.df$Timepoint_ordered) %>% anova()
betadisper(braycurtis_subset, Recipient_data.df$Group) %>% anova()

#Plotting dispersion of Timepoint_ordered
bd_timepoint_ordered <-betadisper(braycurtis_subset, Recipient_data.df$Timepoint_ordered)

boxplot(bd_timepoint_ordered)
plot(bd_timepoint_ordered)

#plotting dispersion of Group
bd_group <- betadisper(braycurtis_subset, Recipient_data.df$Group)

boxplot(bd_group)
plot(bd_group)

#----Bray Curtis pairwise statistics for each timepoint----

#make OTU_raw as a matrix for usability
otu_mat <- as(OTU_raw, "matrix")

#transpose as taxa = rows and samples = columns. Samples should = rows
if(taxa_are_rows(OTU_raw)) {
  otu_mat <- t(otu_mat)
}

#check names match - should return TRUE
all(Recipient_data.df$Sample_ID %in% rownames(otu_mat))

## Stats: (Question: when do groups differ?)
results_list <- list()

for(tp in unique(Recipient_data.df$Timepoint_ordered)) {
  
  subset_df <- Recipient_data.df %>%
    filter(Timepoint_ordered == tp)
  
  if(nrow(subset_df) < 2) next
  if(length(unique(subset_df$Group)) < 2) next
  
  # keep only matching samples
  subset_df <- subset_df[subset_df$Sample_ID %in% rownames(otu_mat), ]
  
  otu_sub <- otu_mat[subset_df$Sample_ID, , drop = FALSE]
  
  dist_sub <- vegdist(otu_sub, method = "bray")
  
  pw <- pairwise.adonis2(
    dist_sub ~ Group,
    data = subset_df,
    permutations = 999
  )
  
  results_list[[tp]] <- pw
}

#bind results for ease of viewing
library(tibble)
pairwise_results <- bind_rows(results_list, .id = "Timepoint")
#view object to see results

#----BC_timepoint_plot - all samples on one graph ----

BC_timepoint_plot <- ggplot(PC_data.df, aes(x = BC1, y = BC2, 
                                  colour = Timepoint_ordered,
                                  size = Group)) + 
  geom_point(alpha = 0.8) +
  scale_size_manual(values = c(
    "AutoFMT"             = 1,
    "Donor match FMT"     = 1,
    "Saline"              = 1,
    "C57/BL6 microbiome"  = 7,
    "LPJ microbiome"      = 7
  ), guide = "none") +  # hide size legend if not needed
  geom_point(
    data = PC_data.df %>% 
      filter(Group %in% c("AutoFMT", "Donor match FMT", "Saline")),
    aes(shape = Group),
    size = 2,          # size has to be +1 from previous to overlay the shape
    stroke = 1
  )+
  scale_color_viridis_d(option = "plasma") +  # discrete palette for ordered factors such as (Timepoint_ordered)
  theme_minimal(base_size = 14) +
  labs(
    title = "Bray–Curtis PCoA (coloured by Timepoint)",
    x = "PCoA Axis 1",
    y = "PCoA Axis 2",
    colour = "Timepoint"
  )

BC_timepoint_plot

#GvHD score mapping to bray curtis
BC_GvHDscore_plot <- ggplot(PC_data.df, aes(x = BC1, y = BC2, 
                                            colour = GvHD_score, #Timepoint_ordered
                                            size = Group)) +
  geom_point(alpha = 0.8) +
  scale_size_manual(values = c(
    "AutoFMT"             = 1,
    "Donor match FMT"     = 1,
    "Saline"              = 1,
    "C57/BL6 microbiome"  = 7,
    "LPJ microbiome"      = 7
  ), guide = "none") +  # hide size legend if not needed
  geom_point(
#    data = PC_data.df %>% 
#      filter(Group %in% c("AutoFMT", "Donor match FMT", "Saline")),
    aes(shape = Group),
    size = 2,          # size has to be +1 from previous to overlay the shape
    stroke = 1
  )+
  scale_color_viridis_c(option = "plasma",
                        breaks = c(0, 1, 2, 3, 4, 5, 6, 7)) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Bray–Curtis PCoA by GvHD-score",
    x = "PCoA Axis 1",
    y = "PCoA Axis 2",
  )

BC_GvHDscore_plot


#graph by study ID

BC_study_plot <- ggplot(PC_data.df, aes(x = BC1, y = BC2, 
                                            colour = Study_ID,
                                            size = Group)) + 
  geom_point(alpha = 0.8) +
  scale_size_manual(values = c(
    "AutoFMT"             = 1,
    "Donor match FMT"     = 1,
    "Saline"              = 1,
    "C57/BL6 microbiome"  = 7,
    "LPJ microbiome"      = 7
  ), guide = "none") +  # hide size legend if not needed
  geom_point(
    data = PC_data.df %>% 
      filter(Group %in% c("AutoFMT", "Donor match FMT", "Saline")),
    aes(shape = Group),
    size = 2,          # size has to be +1 from previous to overlay the shape
    stroke = 1
  )+
  scale_color_viridis_d(option = "plasma") +  # discrete palette for ordered factors
  theme_minimal(base_size = 14) +
  labs(
    title = "Bray–Curtis PCoA (coloured by studyID)",
    x = "PCoA Axis 1",
    y = "PCoA Axis 2",
    colour = "Study ID"
  )

BC_study_plot

#----BC_PLOTS - graphs for each timepoint----
#identify donor samples to overlay over graphs
donor_samples <- c("25-03302", "25-03303")
# Keep only real timepoints (for plotting individual graphs) in order
timepoints <- setdiff(levels(PC_data.df$Timepoint_ordered),
                      c("LPJ microbiome", "C57/BL6 microbiome"))

BC_PLOTS <- list()

for(tp in timepoints){
  df_tp <- PC_data.df %>% filter(Timepoint_ordered == tp)
  
  BC_PLOTS[[tp]] <- ggplot(df_tp, aes(x = BC1, y = BC2, colour = Group, size = Group)) +
    geom_point(aes(shape = Death_timepoint), alpha = 0.8, stroke = 0.5) +
    scale_shape_manual(values = c(
      "endpoint" = 16,            # circle
      "non-endpoint" = 17,        # triangle
      "died preF7" = 15           # square
    )) +
    # ellipses for each group at that timepoint
    stat_ellipse(level = 0.95, size = 0.5) +
    # overlay donor samples: keep group color, but larger + black outline
    geom_point(data = PC_data.df %>% filter(Sample_ID %in% donor_samples),
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

#----Graphing distance to donor microbiome----
#need to assign which sample is which and then generate calculations from those distances and slot into a new column in Recipient.data.df
#plot distance to AutoFMT, DonorFMT and distance to baseline

#create a distance matrix with all samples (includes donors)
bc_full <- as.matrix(braycurtis_dist) 
#technically I have this as braycurtis_matrix already...can clean up this step. 

#define samples to extract distances from
donor_match_id <- "25-03303"
auto_fmt_id    <- "25-03302"

#add calculations of distance to recipient_data.df 
#distance to LPJ donor
Recipient_data.df$dist_to_donor <- bc_full[
  Recipient_data.df$Sample_ID, donor_match_id
]
#distance to C57 donor
Recipient_data.df$dist_to_auto <- bc_full[
  Recipient_data.df$Sample_ID, auto_fmt_id
]

##plot: distance to LPJ donor

distance_to_LPJ <- ggplot(Recipient_data.df, aes(Timepoint_ordered, dist_to_donor,
                              color = Group, group = Mouse_ID)) +
  geom_point(alpha = 0.6) +
  geom_line(alpha = 0.4) +
  stat_summary(aes(group = Group), fun = mean, geom = "line", linewidth = 1.5) +
  theme_minimal() +
  labs (title = "Distance to Donor-match FMT over time",
        y = "Bray-Curtis distance")
distance_to_LPJ

##plot: distance to C57 donor

distance_to_C57 <- ggplot(Recipient_data.df, aes(Timepoint_ordered, dist_to_auto,
                                                 color = Group, group = Mouse_ID)) +
  geom_point(alpha = 0.6) +
  geom_line(alpha = 0.4) +
  stat_summary(aes(group = Group), fun = mean, geom = "line", linewidth = 1.5) +
  theme_minimal() +
  labs (title = "Distance to Auto FMT over time",
        y = "Bray-Curtis distance")
distance_to_C57

#stats
library(lme4)

distance_to_LPJ_model <-lmer(dist_to_donor ~ Group * Timepoint + Study_ID + (1 | Mouse_ID),
                             data = Recipient_data.df)
summary(distance_to_LPJ_model)


distance_to_C57_model <-lmer(dist_to_auto ~ Group * Timepoint + Study_ID + (1 | Mouse_ID),
                             data = Recipient_data.df)
summary(distance_to_C57_model)

#----Graphing distance to baseline microbiome----

#generate distance to baseline calculations

Recipient_data.df$dist_to_baseline <- sapply(Recipient_data.df$Sample_ID, function(sid){ #using sample_ID column instead of rownames to avoid bugs
  
  mouse_id <- Recipient_data.df[sid, "Mouse_ID"]
  
  #find that mouse's baseline sample (S0)
  baseline_sample <-rownames(Recipient_data.df)[
    Recipient_data.df$Mouse_ID == mouse_id &
      Recipient_data.df$Timepoint=="S0"
  ]
  
  #safety check
  
  if (length(baseline_sample) !=1) return(NA)
  
  braycurtis_matrix[sid, baseline_sample]
})

#plot: distance to baseline

distance_to_baseline_plot <-ggplot(Recipient_data.df, aes(Timepoint_ordered, dist_to_baseline, 
                                                            color = Group, group = Mouse_ID)) +
  geom_point(alpha = 0.6) +
  geom_line(alpha = 0.4) +
  stat_summary(aes(group = Group), fun = mean, geom = "line", linewidth = 1.5) +
  theme_minimal()+
  labs (title = "Distance to baseline over time",
        y = "Bray-Curtis distance")

distance_to_baseline_plot

###Stats
distance_to_baseline_model <- lmer(dist_to_baseline ~ Group * Timepoint + Study_ID + (1 | Mouse_ID), 
                                   data = Recipient_data.df)

summary(distance_to_baseline_model)


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

#save plots in order of production
save_plot(Shannon_timeline_plot, "Shannon_timeline_plot_all-mice.pdf")
save_plot(raincloud_shannon_plot, "Raincloud_shannon_plot.pdf")
save_plot(Shannon_all_plot, "Shannon_all_plot.pdf", width = 10, height = 14)
save_plot(BC_timepoint_plot, "BC_timepoint_plot.pdf")
save_plot(BC_GvHDscore_plot, "BC_GvHDscore_plot.pdf")
save_plot(BC_time_plots,"BrayCurtis_PCoA_allTimepoints.pdf", width = 12, height = 15)
save_plot(distance_to_LPJ, "BrayCurtis_distance_to_LPJ.pdf")
save_plot(distance_to_C57, "BrayCurtis_distance_to_C57.pdf")
save_plot(distance_to_baseline_plot, "BrayCurtis_distance_to_baseline.pdf")

