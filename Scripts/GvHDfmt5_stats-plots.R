### Description: Statistics and Graphing 16S dataset      ###
### Author: Anna Li, date: 17092025 - 27082026            ###
### Email: anna.li@adelaide.edu.au                        ###

## Load Libraries
library(ape)              # Required for read.tree
library(phyloseq)         # Required for convert relative abundance to phyloseq object
library(vegan)            # Required for adonis2 - beta diversity
library(ggplot2)          # Required for graphing
library(ggpubr)           # Required for stat_compare_means when graphing shannon
library(gridExtra)        # Required for arranging plots (grid.arrange)
library(dplyr)
library(tidyverse)
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

#Aitchison
ps_clr <- microbiome::transform(physeq_prop, "clr") #input is proportional
clr_counts <- data.frame(
  otu_table(ps_clr),
  check.names = FALSE
)
aitchison_dist = dist(t(clr_counts),method = "euclidean") #The Aitchison distance is just Euclidean on CLR transformed

# MDS (PCoA on Aitchison distance)
mds <- cmdscale(aitchison_dist, k = 2, eig = TRUE)

mds_points <- mds$points
eigvals <- mds$eig

# % variation explained (Aitchison should be Euclidean -> typically all +ve,
# but we guard just in case)
pos_eig <- eigvals[eigvals > 0]
pct_var <- pos_eig / sum(pos_eig) * 100

#define labels
xlab_mds1 <- paste0("MDS1 (", round(pct_var[1], 1), "%)")
ylab_mds2 <- paste0("MDS2 (", round(pct_var[2], 1), "%)")

#----Create a single dataframe/table to reference for graphing----
PC_data.df = data.frame(map_data,alpha_div,
                        "BC1" = braycurtis_pc$vectors[,1],
                        "BC2" = braycurtis_pc$vectors[,2],
                        "UF1" = unifrac_pc$vectors[,1],
                        "UF2" = unifrac_pc$vectors[,2],
                        "WU1" = wunifrac_pc$vectors[,1],
                        "WU2" = wunifrac_pc$vectors[,2],
                        "MDS1" = mds_points[,1],
                        "MDS2" = mds_points[,2])

#Okabe–Ito-derived palette
group_cols <- c(
  "Syn FMT"             = "#D55E00",  # vermillion
  "Donor match FMT"     = "#009E73",  # bluish green
  "Saline"              = "#0072B2",  # blue
  "C57BL/6 microbiome"  = "#E69F00",  # orange 
  "LPJ microbiome"      = "#CC79A7"   # purple-pink 
)

## Adding a "timepoint_ordered" column for longitudinal graphs
# Define the order of the timepoints
timepoint_levels <- c("S0", "D0", "D15", "F7", "F14", "F21", "endpoint", "cull")

# Add an ordered factor column 
PC_data.df$Timepoint_ordered <- factor(PC_data.df$Timepoint,
                                       levels = timepoint_levels,
                                       ordered = TRUE)

#Define C57BL/6 and LPJ microbiome in 'timepoint_ordered'
PC_data.df$Timepoint_ordered <- as.character(PC_data.df$Timepoint_ordered)
PC_data.df$Timepoint_ordered[PC_data.df$Sample_ID %in% c("25-03302")] <- "C57BL/6 microbiome"
PC_data.df$Timepoint_ordered[PC_data.df$Sample_ID %in% c("25-03303")] <- "LPJ microbiome"
PC_data.df$Timepoint_ordered <- factor(PC_data.df$Timepoint_ordered, 
                                       levels = 
                                         c("C57BL/6 microbiome","S0","D0","D15","F7","F14","F21","endpoint","cull","LPJ microbiome"))

#=====Graphing alpha diversity====
# Create object without FMT samples for alpha diversity plots
Recipient_data.df <- PC_data.df[!rownames(PC_data.df)%in% c("25-03302", "25-03303"), ]
# Remove "cull" samples (only plotting stool samples)
Recipient_data.df <- Recipient_data.df %>%
  filter(Timepoint_ordered != "cull")

# Create object where data is only included if Death_timepoint column = endpoint (mouse survived through to end of experiment) and remove cull
Recipient_data_endpoint.df <- Recipient_data.df %>%
  filter(Death_timepoint == "endpoint", Timepoint_ordered != "cull")

#Calculating mean for groups across timepoints 
shannon_summary <- Recipient_data.df %>%
  group_by(Group, Timepoint_ordered) %>%
  summarise(
    mean_shannon = mean(Shannon, na.rm = TRUE),
    se_shannon   = sd(Shannon, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

#Create matrix of timepoint comparisons - for overview of all stats comparisons, select for graphing
timepoint_comparisons <- combn(
  levels(Recipient_data.df$Timepoint_ordered),
  2,
  simplify = FALSE
)

#----Graphing shannon diversity across timepoints----
Shannon_timeline_plot <- ggplot(Recipient_data.df, aes(x = Timepoint_ordered, y = Shannon)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
# geom_line(aes(group = Mouse_ID), alpha = 0.2, linewidth = 0.2) +
# geom_point(aes(shape = Death_timepoint, colour = Group), size = 2, alpha = 0.9) +
  geom_jitter(aes(shape = Death_timepoint, colour = Group), width = 0.2, size = 2, alpha = 0.7) +
  # Mean Shannon line per group
  geom_line(data = shannon_summary,
            aes(x = Timepoint_ordered, y = mean_shannon, group = Group, color = Group),
            linewidth = 1.2) +
  geom_point(data = shannon_summary,
             aes(x = Timepoint_ordered, y = mean_shannon, color = Group),
             size = 3) +
  stat_compare_means(
    comparisons = list(
      c("S0", "D0"),
      c("S0", "D15"),
      c("D0", "D15"),
      c("D15", "F7"),
      c("F7", "F14")
    ),
    method = "wilcox.test",
    p.adjust.method = "holm",
    label = "p.signif", #p.format for specific value
  )+
  scale_color_manual(values = group_cols)+
  scale_shape_manual(values = c(
    "endpoint" = 16,            # circle
    "non-endpoint" = 17,        # triangle
    "died preF7" = 15           # square
  )) +
  theme_minimal() +
  labs(title = "Shannon's diversity per timepoint",
       x = "Timepoint",
       y = "Shannon Diversity Index") +
  theme(legend.position = "right")

Shannon_timeline_plot

#----Graphing Shannon's diversity as a raincloud----
library(ggdist)
library(gghalves)

#Plot: raincloud_shannon
raincloud_shannon_plot <- ggplot(Recipient_data.df, 
                                 aes(x = Timepoint_ordered, y = Shannon)) +
  # Half-violin (density)
  stat_halfeye(
    adjust = 0.7, width = 0.7, .width = 0, justification = -0.2,
    point_interval = median_qi,
    alpha = 0.6) +
  # Jittered points
  geom_point(
    aes(shape = Death_timepoint, color = Group), 
    position = position_jitter(width = 0.2, height = 0),
    alpha = 0.6, size = 2.5) +
  # Boxplot (optional; remove if you want summary dot only)
  geom_boxplot(
    width = 0.1, outlier.shape = NA,
    alpha = 0.7, color = "black") +
  # Mean Shannon line per group
  geom_line(data = shannon_summary,
            aes(x = Timepoint_ordered, y = mean_shannon, group = Group, color = Group),
            linewidth = 1.2) +
  geom_point(data = shannon_summary,
             aes(x = Timepoint_ordered, y = mean_shannon, color = Group),
             size = 3) +
  stat_compare_means(
    comparisons = list(
      c("S0", "D0"),
      c("S0", "D15"),
      c("D0", "D15"),
      c("D15", "F7"),
      c("F7", "F14")
    ),
    method = "wilcox.test",
    p.adjust.method = "holm",
    label = "p.signif", #p.format for specific value
  )+
  scale_color_manual(values = group_cols)+
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 14,
                              face = "bold",
                              hjust = 0.5)) +
  labs(title = "Shannon Diversity by Timepoint",
       x = "",
       y = "Shannon Diversity Index")

raincloud_shannon_plot

#----Graphing Shannon's diversity - by group at specific timepoints----
# Generate all pairwise combinations of groups/define pairwise comparisons
comparisons <- combn(unique(Recipient_data.df$Group), 2, simplify = FALSE)

#Plotting Shannon's diversity - at D15 by group
# Subset only D15 samples
alpha_div_D15 <- Recipient_data.df %>%
  filter(Timepoint_ordered == "D15")
# Kruskal-Wallis test across all groups at endpoint
kruskal_result <- kruskal.test(Shannon ~ Group, data = alpha_div_D15)
print(kruskal_result)

# Boxplot with kruskal-Wallis comparisons and wilcox
shannon_D15_plot <- ggplot(alpha_div_D15, aes(x = Group, y = Shannon, fill = Group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA,) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7) +
  scale_fill_manual(values = group_cols)+
  theme_minimal(base_size = 14) +
  labs(title = "Shannon Diversity at D15 by Group",
       x = "",
       y = "Shannon Diversity Index") +
  guides(fill = "none") +
  theme(legend.position = "none") +
  stat_compare_means(method = "kruskal.test", label.y = max(alpha_div_D15$Shannon) + 0.2) + # global test
  stat_compare_means(method = "wilcox.test", 
                     comparisons = combn(unique(alpha_div_D15$Group), 2, simplify = FALSE),
                    label = "p.signif")
shannon_D15_plot

#----Graphing Shannon's diversity - by group for all individual timepoints---- 

#stats
library(rstatix)

wilcox_results <- Recipient_data.df %>%
  group_by(Timepoint_ordered) %>%
  pairwise_wilcox_test(
    Shannon ~ Group,
    p.adjust.method = "BH"
  )

wilcox_results

wilcox_results <- wilcox_results %>%
  add_y_position(
    fun = "max",
    step.increase = 0.2)

#  stat_compare_means(method = "wilcox.test",
comparisons = combn(unique(Recipient_data.df$Group), 2, simplify = FALSE),
label = "p.signif") +

Shannon_all_plot <- ggplot(Recipient_data.df, aes(x = Group, y = Shannon, fill = Group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  scale_fill_manual(values = group_cols) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7) +
  stat_pvalue_manual(wilcox_results,
                     label = "p.adj.signif",
                     hide.ns = TRUE) +
  facet_wrap(~ Timepoint_ordered, nrow=2) +   # cycle through timepoints
  theme_classic2(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 0),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 12,
                              face = "bold",
                              hjust = 0.5),
    axis.line = element_line(
      linewidth = 0.5,
      linetype = "solid"
    )) +
  labs(title = "Shannon's diversity across timepoints by treatment group",
       x = "",
       y = "Shannon Diversity Index")
#  guides(fill = "none") +
  theme(legend.position = "inside")
Shannon_all_plot

#demonstrates that there are no statistic differences in alpha diversity at each timepoint

#====Graphing beta diversity====
#----Bray Curtis and Aitchison global statistics / PERMANOVA----
#Note: use recipient_data for stats (don't factor in donor samples for global stats on trajectory yet)
#Align Recipient_data.df with braycurtis_distance 
braycurtis_matrix <- as.matrix(braycurtis_dist)
samples <- intersect(rownames(braycurtis_matrix), rownames(Recipient_data.df))

Recipient_data.df <- Recipient_data.df[samples, ]
braycurtis_subset <- as.dist(braycurtis_matrix[samples, samples])

#Repeat for Aitchison - recipient_data.df with aitchison distance 
aitchison_matrix<- as.matrix(aitchison_dist)
samples <- intersect(rownames(aitchison_matrix), rownames(Recipient_data.df))

Recipient_data.df <- Recipient_data.df[samples, ]
aitchison_subset <- as.dist(aitchison_matrix[samples, samples])

# Make sure that timepoint_ordered is a factor not a numeric ^should have been coded before but redo here
Recipient_data.df$Timepoint_ordered <- factor(
  Recipient_data.df$Timepoint_ordered,
  levels = c("S0","D0","D15","F7","F14","F21","endpoint"),
  ordered = TRUE
)

###Stats: (Question: is variation in gut microbiome explained by time?/Do microbiome compositions differ between timepoints?)
adonis2(aitchison_subset ~ Timepoint_ordered, 
        data = Recipient_data.df,
        permutations = 999,
        strata = Recipient_data.df$Mouse_ID) #permutations are constrained within each mouse/repeated measures

###Stats: (Question: does synFMT vs Donor-matchFMT vs Saline lead to different longitudinal changes?)
adonis2(aitchison_subset ~ Group * Timepoint_ordered,
        data = Recipient_data.df,
        permutations = 999,
        strata = Recipient_data.df$Mouse_ID,
    by = "margin") #term level output

#checking inter dispersion within the time and group (PERMANOVA detects group centroids i.e. biological or dispersion i.e. variances)
betadisper(aitchison_subset, Recipient_data.df$Timepoint_ordered) %>% anova()
betadisper(aitchison_subset, Recipient_data.df$Group) %>% anova()

#Plotting dispersion of Timepoint_ordered
bd_timepoint_ordered <-betadisper(braycurtis_subset, Recipient_data.df$Timepoint_ordered)
aitch_timepoint_ordered <-betadisper(aitchison_subset, Recipient_data.df$Timepoint_ordered)

boxplot(bd_timepoint_ordered)
plot(bd_timepoint_ordered)

boxplot(aitch_timepoint_ordered)
plot(aitch_timepoint_ordered)

#plotting dispersion of Group
bd_group <- betadisper(braycurtis_subset, Recipient_data.df$Group)
aitch_group <- betadisper(aitchison_subset, Recipient_data.df$Group)

boxplot(bd_group)
plot(bd_group)

boxplot(aitch_group)
plot(aitch_group)

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

#----Aitchison pairwise statistics for each timepoint----

#stats: (Question: how do groups differ?)

# Convert aitchison_subset (contains only recipient samples) into a matrix for usability
aitchison_subset_mat<-as.matrix(aitchison_subset)
#check names match - should return TRUE
all(Recipient_data.df$Sample_ID %in% rownames(aitchison_subset_mat))
all(rownames(Recipient_data.df)%in% rownames(aitchison_subset_mat))

results_list_a <- list()

for(tp in unique(Recipient_data.df$Timepoint_ordered)) {
  
  subset_df <- Recipient_data.df %>%
    filter(Timepoint_ordered == tp)
  
  if(nrow(subset_df) < 2) next
  if(length(unique(subset_df$Group)) < 2) next
  
  sample_ids <- subset_df$Sample_ID
  
  dist_sub <- as.dist(
    aitchison_subset_mat[sample_ids, sample_ids]
  )
  
  pw <- pairwise.adonis2(
    dist_sub ~ Group,
    data = subset_df,
    permutations = 999
  )
  
  results_list_a[[as.character(tp)]] <- pw
}

pairwise_results_a <- bind_rows(results_list, .id = "Timepoint")


#----BC_timepoint_plot - all samples on one graph ----

BC_timepoint_plot <- ggplot(PC_data.df, aes(x = BC1, y = BC2, 
                                  colour = Timepoint_ordered,
                                  size = Group)) + 
  geom_point(alpha = 0.8) +
  scale_size_manual(values = c(
    "Syn FMT"             = 1,
    "Donor match FMT"     = 1,
    "Saline"              = 1,
    "C57BL/6 microbiome"  = 7,
    "LPJ microbiome"      = 7
  ), guide = "none") +  # hide size legend if not needed
  geom_point(
    data = PC_data.df %>% 
      filter(Group %in% c("Syn FMT","Donor match FMT","Saline")),
    aes(shape = Group),
    size = 2,          # size has to be +1 from previous to overlay the shape
    stroke = 1
  )+
  scale_shape_manual(values = c(
    "Donor match FMT" = 16,            # circle
    "Saline" = 17,                     # triangle
    "Syn FMT" = 15                     # square
  )) +
  scale_color_viridis_d(option = "plasma", direction = -1) +  # discrete palette for ordered factors such as (Timepoint_ordered)
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
    "Syn FMT"             = 1,
    "Donor match FMT"     = 1,
    "Saline"              = 1,
    "C57BL/6 microbiome"  = 7,
    "LPJ microbiome"      = 7
  ), guide = "none") +  # hide size legend if not needed
  geom_point(
#    data = PC_data.df %>% 
#      filter(Group %in% c("Syn FMT", "Donor match FMT", "Saline")),
    aes(shape = Group),
    size = 2,          # size has to be +1 from previous to overlay the shape
    stroke = 1
  )+
  scale_shape_manual(values = c(
    "Donor match FMT" = 16,            # circle
    "Saline" = 17,                     # triangle
    "Syn FMT" = 15,                     # square
    "LPJ microbiome" = 1,
    "C57BL/6 microbiome" = 0
  )) +
  scale_color_viridis_c(option = "plasma", direction = -1,
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
    "Syn FMT"             = 1,
    "Donor match FMT"     = 1,
    "Saline"              = 1,
    "C57BL/6 microbiome"  = 7,
    "LPJ microbiome"      = 7
  ), guide = "none") +  # hide size legend if not needed
  geom_point(
    data = PC_data.df %>% 
      filter(Group %in% c("Syn FMT", "Donor match FMT", "Saline")),
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
                      c("LPJ microbiome", "C57BL/6 microbiome"))

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
      "Syn FMT"            = 2,
      "Donor match FMT"    = 2,
      "Saline"             = 2,
      "C57BL/6 microbiome" = 6,
      "LPJ microbiome"     = 6
    ), guide = "none") +
    theme_minimal(base_size = 14) +
    labs(title = paste("Bray–Curtis PCoA:", tp),
         x = "PCoA Axis 1",
         y = "PCoA Axis 2")
}


library(cowplot)
BC_time_plots <- plot_grid(plotlist = BC_PLOTS[timepoints],
                           ncol = 2,
                           labels = c("i", "ii", "iii", "iv", "v", "vi", "Vii","viii"))
BC_time_plots

#----Aitchison MDS plots----

MDS_timepoint_plot <- ggplot(PC_data.df, aes(x = MDS1, y = MDS2, 
                                             colour = Timepoint_ordered,
                                             size = Group)) + 
  geom_point(alpha = 0.8) +
  scale_size_manual(values = c(
    "Syn FMT"             = 1,
    "Donor match FMT"     = 1,
    "Saline"              = 1,
    "C57BL/6 microbiome"  = 7,
    "LPJ microbiome"      = 7
  ), guide = "none") +  # hide size legend if not needed
  geom_point(
    data = PC_data.df %>% 
      filter(Group %in% c("Syn FMT","Donor match FMT","Saline")),
    aes(shape = Group),
    size = 2,          # size has to be +1 from previous to overlay the shape
    stroke = 1
  )+
  scale_shape_manual(values = c(
    "Donor match FMT" = 16,            # circle
    "Saline" = 17,                     # triangle
    "Syn FMT" = 15                     # square
  )) +
  scale_color_viridis_d(option = "plasma", direction = -1) +  # discrete palette for ordered factors such as (Timepoint_ordered)
  theme_minimal(base_size = 14) +
  labs(
    title = "Aitchisons MDS (coloured by Timepoint)",
    x = xlab_mds1,
    y = ylab_mds2,
    colour = "Timepoint"
  )

MDS_timepoint_plot

#GvHD score mapping to MDS

MDS_GvHDscore_plot <- ggplot(PC_data.df, aes(x = MDS1, y = MDS2, 
                                            colour = GvHD_score, #Timepoint_ordered
                                            size = Group)) +
  geom_point(alpha = 0.8) +
  scale_size_manual(values = c(
    "Syn FMT"             = 1,
    "Donor match FMT"     = 1,
    "Saline"              = 1,
    "C57BL/6 microbiome"  = 7,
    "LPJ microbiome"      = 7
  ), guide = "none") +  # hide size legend if not needed
  geom_point(
    #    data = PC_data.df %>% 
    #      filter(Group %in% c("Syn FMT", "Donor match FMT", "Saline")),
    aes(shape = Group),
    size = 2,          # size has to be +1 from previous to overlay the shape
    stroke = 1
  )+
  scale_shape_manual(values = c(
    "Donor match FMT" = 16,            # circle
    "Saline" = 17,                     # triangle
    "Syn FMT" = 15,                     # square
    "LPJ microbiome" = 1,
    "C57BL/6 microbiome" = 0
  )) +
  scale_color_viridis_c(option = "plasma", direction = -1,
                        breaks = c(0, 1, 2, 3, 4, 5, 6, 7)) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Aitchisions MDS by GvHD-score",
    x = xlab_mds1,
    y = ylab_mds2
  )

MDS_GvHDscore_plot

#MDS plots for each timepoint

MDS_PLOTS <- list()

for(tp in timepoints){
  df_tp <- PC_data.df %>% filter(Timepoint_ordered == tp)
  
  MDS_PLOTS[[tp]] <- ggplot(df_tp, aes(x = MDS1, y = MDS2, colour = Group, size = Group)) +
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
               aes(x = MDS1, y = MDS2, colour = Group), 
               size = 5, shape = 21, stroke = 1.5, fill = NA) +
    scale_color_manual(values = group_cols) +
    scale_size_manual(values = c(
      "Syn FMT"            = 2,
      "Donor match FMT"    = 2,
      "Saline"             = 2,
      "C57BL/6 microbiome" = 6,
      "LPJ microbiome"     = 6
    ), guide = "none") +
    theme_bw(base_size = 14) +
    theme(legend.position="none")+
    labs(title = paste(tp),
         x = xlab_mds1,
         y = ylab_mds2)
}

legend_plot <- MDS_PLOTS[[1]] +
  theme(
    legend.position = "bottom",
  )
legend_plot
legend <- get_legend(legend_plot)

library(cowplot)
MDS_time_grid <- plot_grid(plotlist = MDS_PLOTS[timepoints],
                           ncol = 4,
                           labels = c("i", "ii", "iii", "iv", "v", "vi", "Vii","viii"))
MDS_time_plots <- plot_grid(
  MDS_time_grid,
  legend,
  ncol = 1,
  rel_heights = c(1, 0.08)
)
MDS_time_plots

#----Graphing distance to donor microbiome----
#need to assign which sample is which and then generate calculations from those distances and slot into a new column in Recipient.data.df
#plot distance to Syn FMT, Donor match FMT and distance to baseline

#create distance matrix with all samples (includes donors)
braycurtis_matrix <- as.matrix(braycurtis_dist) #already completed in above script but run again
aitchison_matrix<- as.matrix(aitchison_dist) #already completed in above script but run again

#define samples to extract distances from
donor_match_id <- "25-03303"
syn_fmt_id    <- "25-03302"

#add calculations of distance to recipient_data.df 
#distance to LPJ donor
Recipient_data.df$dist_to_donor <- braycurtis_matrix[
  Recipient_data.df$Sample_ID, donor_match_id
]

Recipient_data.df$dist_to_donor_a <- aitchison_matrix[
  Recipient_data.df$Sample_ID, donor_match_id
]

#distance to C57 donor
Recipient_data.df$dist_to_syn <- braycurtis_matrix[
  Recipient_data.df$Sample_ID, syn_fmt_id
]

Recipient_data.df$dist_to_syn_a <- aitchison_matrix[
  Recipient_data.df$Sample_ID, syn_fmt_id
]

##plot: distance to LPJ donor

distance_to_LPJ <- ggplot(Recipient_data.df, aes(Timepoint_ordered, dist_to_donor_a, #swap between dist_to_donor_a***
                              color = Group, group = Mouse_ID)) +
  geom_line(alpha = 0.4) +
  stat_summary(aes(group = Group), fun = mean, geom = "line", linewidth = 1.5) +
  scale_color_manual(values = group_cols) +
  theme_minimal(base_size =16) +
  theme(
    plot.title = element_text(size = 16,
                              face = "bold",
                              hjust = 0.5),
    axis.line = element_line(linewidth = 0.5)
  ) +
  labs (title = "Distance to donor match (LPJ) FMT over time",
        y = "Aitchison distance",
        x = "")
distance_to_LPJ

##plot: distance to C57 donor

distance_to_C57 <- ggplot(Recipient_data.df, aes(Timepoint_ordered, dist_to_syn_a,
                                                 color = Group, group = Mouse_ID)) +
  geom_line(alpha = 0.4) +
  stat_summary(aes(group = Group), fun = mean, geom = "line", linewidth = 1.5) +
  scale_color_manual(values = group_cols) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(size = 16,
                              face = "bold",
                              hjust = 0.5),
    axis.line = element_line(linewidth = 0.5)
  ) +
  labs (title = "Distance to Syn (C57BL/6) FMT over time",
        y = "Aitchison distance",
        x = "")
distance_to_C57

###stats

#Create mixed model of distance to LPJ 
library(lme4)
distance_to_LPJ_model <-lmer(dist_to_donor_a ~ Group * Timepoint + Study_ID + (1 | Mouse_ID),
                             data = Recipient_data.df)
summary(distance_to_LPJ_model)

#Overall effects of group x timepoint based on model
anova(distance_to_LPJ_model)

#Estimate group differences at each timepoint based on model
library(emmeans)
emmeans(distance_to_LPJ_model, 
        pairwise ~ Group | Timepoint,
        adjust = "BH")

#Create mixed model of distance to C57
distance_to_C57_model <-lmer(dist_to_syn_a ~ Group * Timepoint + Study_ID + (1 | Mouse_ID),
                             data = Recipient_data.df)
summary(distance_to_C57_model)

#Overall effects of group x timepoint based on model
anova(distance_to_C57_model)

#Estimate group differences at each timepoint based on model
emmeans(distance_to_C57_model, 
        pairwise ~ Group | Timepoint,
        adjust = "BH")

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
  
  aitchison_matrix[sid, baseline_sample] ##switch out between matrix of bray or aitchison here***
})

#plot: distance to baseline

distance_to_baseline_plot <-ggplot(Recipient_data.df, aes(Timepoint_ordered, dist_to_baseline, 
                                                            color = Group, group = Mouse_ID)) +
  geom_line(alpha = 0.4) +
  stat_summary(aes(group = Group), fun = mean, geom = "line", linewidth = 1.5) +
  scale_color_manual(values = group_cols) +
  theme_minimal(base_size = 16)+
  theme(
    plot.title = element_text(size = 16,
                              face = "bold",
                              hjust = 0.5),
    axis.line = element_line(linewidth = 0.5)
  ) +
  labs (title = "Distance to baseline over time",
        y = "Aitchison distance",
        x = "")

distance_to_baseline_plot

###Stats
distance_to_baseline_model <- lmer(dist_to_baseline ~ Group * Timepoint + Study_ID + (1 | Mouse_ID), 
                                   data = Recipient_data.df)

#summary(distance_to_baseline_model)
anova(distance_to_baseline_model)

emmeans(
  distance_to_baseline_model,
  pairwise ~ Group | Timepoint,
  adjust = "BH"
)

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
save_plot(Shannon_timeline_plot, "Shannon_timeline_plot_all-mice.pdf", width = 8, height =5)
save_plot(raincloud_shannon_plot, "Raincloud_shannon_plot.pdf", width = 8, height = 5)
save_plot(Shannon_all_plot, "Shannon_all_plot.pdf", width = 15, height = 8)
save_plot(shannon_D15_plot, "shannon_D15_plot.pdf")

save_plot(BC_timepoint_plot, "BC_timepoint_plot.pdf")
save_plot(BC_GvHDscore_plot, "BC_GvHDscore_plot.pdf")
save_plot(BC_time_plots,"BrayCurtis_PCoA_allTimepoints.pdf", width = 12, height = 15)

save_plot(MDS_timepoint_plot, "MDS_timepoint_plot.pdf", width = 7, height = 6)
save_plot(MDS_GvHDscore_plot, "MDS_GvHDscore_plot.pdf", width = 7, height = 6)
save_plot(MDS_time_plots, "Aitchisons_MDS_allTimepoints.pdf", width = 14, height = 8)

save_plot(distance_to_LPJ, "Aitchison_distance_to_LPJ.pdf", width = 8, height= 4)
save_plot(distance_to_C57, "Aitchison_distance_to_C57.pdf", width = 8, height= 4)
save_plot(distance_to_baseline_plot, "Aitchison_distance_to_baseline.pdf", width = 8, height= 4)

