### Description: Graphing taxonomic data                  ###
### Author: Anna Li, date: 04102025                       ###
### Email: anna.li@adelaide.edu.au                        ###

## Load Libraries
library(ape)              # Required for read.tree
library(phyloseq)         # Required for convert relative abundance to phyloseq object
library(vegan)            # Required for adonis2 - beta diversity
library(ggplot2)          # Required for graphing
library(ggpubr)           # Required for stat_compare_means when graphing shannon
library(gridExtra)        # Required for arranging plots (grid.arrange)

#====Load filtered counts, taxonomy, metadata, phyloseq objects====
raw_counts_keep <- readRDS("raw_counts_keep.rds")
taxa_filtered <- readRDS("taxa_filtered.rds")             #will need this
metadata_keep <- readRDS("metadata_keep.rds")
OTU_raw <- readRDS("OTU_raw.rds")
tax_mat<-readRDS("tax_mat.rds")
map_data<-readRDS("map_data.rds")
physeq_raw <- readRDS("physeq_raw.rds")
physeq_prop <- readRDS("physeq_prop.rds")                 #will need this
physeq_raw_tree <-readRDS("physeq_raw_tree.rds")
# plus PC_data.df
PC_data.df <-readRDS("PC_data.df.rds")


##====Creating taxonomic bar plots for donor samples====
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)   # for colour shading functions

# Choose taxonomic level
level <- "Genus"

# Collapse counts to genus
counts_genus <- counts_prop %>%
  mutate(!!level := taxa_filtered[[level]]) %>%     # add level column
  group_by(!!sym(level)) %>%                        # group by level
  summarise(across(everything(), sum)) %>%          # sum across ASVs
  as.data.frame()

#Organise new data object
rownames(counts_genus) <- counts_genus[[level]]
counts_genus[[level]] <- NULL

#plot 
FMT_to_plot <- c("25-03302", "25-03303")
counts_genus_subset <- counts_genus[, FMT_to_plot]

#reshape for ggplot2 to read (long form = one row per observation)
counts_genus_long <- counts_genus_subset %>%
  tibble::rownames_to_column("Taxon") %>%
  pivot_longer(-Taxon, names_to = "SampleID", values_to = "Abundance")

# Example counts_genus_long: columns = Taxon, SampleID, Abundance
# And you also have a taxonomy dataframe with Genus + Class

# Make sure each genus has only ONE class assigned (currently multiple ASV annotated to the same genus)
genus_class <- taxa_filtered %>%
  select(Genus, Class) %>%
  distinct()

# Now join safely
counts_taxa_long <- counts_genus_long %>%
  left_join(genus_class, by = c("Taxon" = "Genus"))

# Build genus-level colour map
genus_list <- counts_taxa_long %>%
  group_by(Class) %>%
  reframe(Genera = unique(Taxon))

library(RColorBrewer)

genus_colors <- unlist(lapply(split(genus_list$Genera, genus_list$Class), function(genera){
  n <- length(genera)
  cols <- colorRampPalette(brewer.pal(min(8, n), "Set3"))(n) # dynamically generate n colours
  names(cols) <- genera
  cols
}))

#checking step
unique(counts_taxa_long$Taxon)
names(genus_colors)
setdiff(unique(counts_taxa_long$Taxon), names(genus_colors))

# FIX: If names are like "Bacilli.Anaeroplasma", keep only the genus
names(genus_colors) <- sapply(strsplit(names(genus_colors), "\\."), `[`, 2)

# plot the graph
Donor_genus_plot <- ggplot(counts_taxa_long, aes(x = SampleID, y = Abundance, fill = Taxon)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = genus_colors) +
  theme_minimal(base_size = 14) +
  labs(title = "Taxonomic composition")

#====Defining unique genera====
library(dplyr)

# Keep only taxa with non-zero abundance in each sample
#autofmt
g1 <- counts_taxa_long %>% 
  filter(SampleID == "25-03302", Abundance > 0) %>% 
  pull(Taxon) %>% 
  unique()
#donorfmt
g2 <- counts_taxa_long %>% 
  filter(SampleID == "25-03303", Abundance > 0) %>% 
  pull(Taxon) %>% 
  unique()

# Genera unique to each sample
unique_g1 <- setdiff(g1, g2)
unique_g2 <- setdiff(g2, g1)

# Genera shared between the two samples
shared_g <- intersect(g1, g2)

# View results
unique_g1
unique_g2
shared_g

#====Defining differentially abundant genera====
library(dplyr)
library(tidyr)

# Define your samples
s1 <- "25-03302" #autoFMT
s2 <- "25-03303" #donorFMT

# Define parameters
fold_cutoff <- 10       # how much more abundant to be considered "distinctive"
min_abundance <- 0.0001 # only consider taxa with at least this abundance
pseudocount <- 1e-6     # avoids divide-by-zero problems

# Convert to wide format
wide_df <- counts_genus_long %>%
  filter(SampleID %in% c(s1, s2)) %>%
  pivot_wider(names_from = SampleID, values_from = Abundance, values_fill = 0)

# Compute fold-change ratios
wide_df <- wide_df %>%
  mutate(
    fold_s1_over_s2 = (!!sym(s1) + pseudocount) / (!!sym(s2) + pseudocount),
    fold_s2_over_s1 = (!!sym(s2) + pseudocount) / (!!sym(s1) + pseudocount)
  )

# Identify "distinctive" taxa for each sample
distinct_s1 <- wide_df %>%
  filter(!!sym(s1) > min_abundance, fold_s1_over_s2 >= fold_cutoff) %>%
  pull(Taxon)

distinct_s2 <- wide_df %>%
  filter(!!sym(s2) > min_abundance, fold_s2_over_s1 >= fold_cutoff) %>%
  pull(Taxon)

# Results
cat("Distinctive taxa for", s1, ":\n")
print(distinct_s1)

cat("\nDistinctive taxa for", s2, ":\n")
print(distinct_s2)

##====Taxonomic plot highlighting unique/differentially abundant genera====
library(dplyr)
library(ggplot2)

# Filter only samples of interest
plot_data <- counts_taxa_long %>%
  filter(SampleID %in% c("25-03302", "25-03303"))

# Calculate label positions per sample
plot_data <- plot_data %>%
  group_by(SampleID) %>%
  arrange(SampleID, Taxon) %>%
  mutate(ypos = cumsum(Abundance) - 0.5 * Abundance)

# Add labels only for unique genera !!OR!!===>
plot_labels <- plot_data %>%
  filter((SampleID == "25-03302" & Taxon %in% unique_g1) |
           (SampleID == "25-03303" & Taxon %in% unique_g2))
# Add labels only for distinct genera
#plot_labels <- plot_data %>%
  #filter((SampleID == "25-03302" & Taxon %in% distinct_s1) |
           #(SampleID == "25-03303" & Taxon %in% distinct_s2))

# changing order to match bar plot and legend 
plot_data$Taxon <- factor(plot_data$Taxon, levels = rev(unique(plot_data$Taxon)))

# Plot with precise label alignment
Donor_genus_plot2 <- ggplot(plot_data, aes(x = SampleID, y = Abundance, fill = Taxon)) +
  geom_bar(stat = "identity", color = "grey30", linewidth = 0.2) +
  geom_text(data = plot_labels, 
            aes(y = ypos, label = Taxon),
            size = 3, color = "black", fontface = "bold") +
  scale_fill_manual(values = genus_colors) +
  theme_minimal(base_size = 14) +
  scale_x_discrete(labels = c("25-03302" = "C57BL/6 microbiome",
                              "25-03303" = "LPJ microbiome")) +
  labs(title = "Taxonomic Composition with Unique Genera Labeled",
       x = "",
       y = "Relative Abundance")

Donor_genus_plot2

#====Plotting changes in unique genera over time for each mouse ====

#Unique_g1
#Unique_g2
#Utilising physeq_raw that contains taxa and mapdata

library(phyloseq)
library(dplyr)
library(tidyr)
library(ggplot2)

# Agglomerate to genus level (or another rank)
physeq_genus <- tax_glom(physeq_raw, taxrank = "Genus", NArm = TRUE)

# Transform to relative abundance
physeq_genus_rel <- transform_sample_counts(physeq_genus, function(x) x / sum(x))

# Convert to a long data frame for plotting
genus_long <- psmelt(physeq_genus_rel)
# genus_long now has:
#   Sample, Genus, Abundance, Mouse_ID, Timepoint_ordered, Group, etc.

# Filter for genera of interest (unique_g1 and unique_g2)
genus_long_unique <- genus_long %>%
  filter(Genus %in% c(unique_g1, unique_g2))

# Add timepoint_ordered
genus_long_unique <- genus_long_unique %>%
  mutate(Timepoint_ordered = factor(Timepoint,
                                    levels = timepoint_levels,
                                    ordered = TRUE))
# filter out unwanted timepoints cecal and NA
genus_long_unique_filt <- genus_long_unique %>%
  filter(!is.na(Timepoint_ordered) & Timepoint_ordered != "cecal")

# Mark which genera belong to unique_g1 and unique_g2
genus_long_unique_filt <- genus_long_unique_filt %>%
  mutate(Origin = case_when(
    Genus %in% unique_g1 ~ "C57BL/6 unique",
    Genus %in% unique_g2 ~ "LPJ unique",
    TRUE ~ "Shared"
  ))

# Plot change over time for all unique genera
All_unique_plot = ggplot(genus_long_unique_filt, aes(x = Timepoint_ordered, y = Abundance,
                              fill = Genus, group = Genus)) +
  geom_area(position = "stack", alpha = 0.9) +
  facet_wrap(~ Group + Mouse_ID, scales = "fixed", ncol = 5) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Temporal dynamics of distinctive genera",
    x = "Timepoint",
    y = "Relative abundance"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#plot LPJ unique genera
LPJ_unique_df <- genus_long_unique_filt %>%
  filter(Genus %in% unique_g2)

LPJ_unique_plot = ggplot(LPJ_unique_df, aes(x = Timepoint_ordered, y = Abundance,
                                   fill = Genus, group = Genus)) +
  geom_area(position = "stack", alpha = 0.9) +
  facet_wrap(~ Group + Mouse_ID, scales = "fixed", ncol = 5) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Temporal dynamics of LPJ unique genera",
    x = "Timepoint",
    y = "Relative abundance"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#plot C57BL/6 unique genera
C57_unique_df <- genus_long_unique_filt %>%
  filter(Genus %in% unique_g1)

C57_unique_plot = ggplot(C57_unique_df, aes(x = Timepoint_ordered, y = Abundance,
                          fill = Genus, group = Genus)) +
  geom_area(position = "stack", alpha = 0.9) +
  facet_wrap(~Group + Mouse_ID, scales = "fixed", ncol = 5) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Temporal dynamics of C57 unique genera",
    x = "Timepoint",
    y = "Relative abundance"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#====Save graphs for export====
# use width 1300 and height 700 save as TIFF
# or directly save from zoom page right click

Donor_genus_plot
Donor_genus_plot2

All_unique_plot
LPJ_unique_plot
C57_unique_plot



