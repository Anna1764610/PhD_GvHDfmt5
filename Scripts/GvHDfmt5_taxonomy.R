### Description: Graphing taxonomy                        ###
### Author: Anna Li, date: 01052026                       ###
### Email: anna.li@adelaide.edu.au                        ###

## Load Libraries
library(ape)              # Required for read.tree
library(phyloseq)         # Required for convert relative abundance to phyloseq object
library(vegan)            # Required for adonis2 - beta diversity
library(ggplot2)          # Required for graphing
library(ggpubr)           # Required for stat_compare_means when graphing shannon
library(gridExtra)        # Required for arranging plots (grid.arrange)
library(dplyr)
library(tidyr)

## Set up folders
project_dir <- here::here()
results_dir <- file.path(project_dir, "Results")

#====Plotting at family level====
#list of families
unique_families <- unique(taxa_filtered$Family)
unique_families

#Assign colours
family_colour <-c(
  "Lactobacillaceae" = "red4",
  "Prevotellaceae" = "blue4",
  "Bacteroidaceae" = "forestgreen",
  "Muribaculaceae" = "goldenrod2",
  "Lachnospiraceae" = "salmon2",
  "Sutterellaceae" = "powderblue",
  "Ruminococcaceae" = "#3C280D",
  "Erysipelatoclostridiaceae" = "deeppink3",
  "Erysipelotrichaceae" = "orchid4",
  "Acholeplasmataceae" = "turquoise4",
  "Clostridia_UCG-014" = "slateblue3",
  "Oscillospiraceae" = "peachpuff1",
  "Enterobacteriaceae" = "lightskyblue4",
  "Staphylococcaceae" = "firebrick3",
  "Anaerovoracaceae" = "springgreen4",
  "Clostridiaceae" = "tan4",
  "Butyricicoccaceae" = "seagreen3",
  "Eggerthellaceae" = "wheat4",
  "Clostridia_vadinBB60_group" = "darkgoldenrod3",
  "Akkermansiaceae" = "dodgerblue3",
  "Atopobiaceae" = "coral4",
  "UCG-010" = "mediumpurple4",
  "Enterococcaceae" = "violetred4", #end of shared family between MBC study and GVHD study
  "[Eubacterium]_coprostanoligenes_group"= "navy",
  "Bifidobacteriaceae" = "aquamarine4",
  "Tannerellaceae"="palevioletred4",
  "Monoglobaceae" =  "indianred3",  
  "Defluviitaleaceae"="khaki4",              
  "Corynebacteriaceae"= "darkseagreen",
  "Gastranaerophilales"= "olivedrab4",                                  
  "Rikenellaceae"="sienna3"
)

#----Create a subset of physeq object physeq_prop with only the donor samples----
donor_samples <- c("25-03303","25-03302")
physeq_donor <- prune_samples(donor_samples, physeq_prop)
#check
sample_names(physeq_donor)

barplot_donor <- tax_glom(physeq_donor, taxrank = "Family")
barplot_donor_family <- psmelt(barplot_donor) #Melt physeq table into a data frame for plotting 

barplot_donor_family$Family[is.na(barplot_donor_family$Family)] <- "Unclassified"



#----Plot donor samples family bar graph----

barplot_donor <- ggplot(barplot_donor_family, aes(x = Mouse_ID, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = family_colour) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(title = paste("Relative Abundance (Family level)"), x = "sample", y = "Relative Abundance")

barplot_donor


#----Create a subset of physeq object physeq_prop with all recipient timepoints====
#define timepoints
Timepoint_ordered = c("S0","D0","D15","F7","F14","F21","Cull")
#fixing physeq_prop phyloseq object to have ordered timepoints
sample_data(physeq_prop)$Timepoint <- factor(
  sample_data(physeq_prop)$Timepoint,
  levels = Timepoint_ordered,
  ordered = TRUE
)

physeq_S0 <- subset_samples(physeq_prop, Timepoint == "S0")
physeq_D0 <- subset_samples(physeq_prop, Timepoint == "D0")
physeq_D15 <- subset_samples(physeq_prop, Timepoint == "D15")
physeq_F7 <- subset_samples(physeq_prop, Timepoint == "F7")
physeq_F14 <- subset_samples(physeq_prop, Timepoint == "F14")
physeq_F21 <- subset_samples(physeq_prop, Timepoint == "F21")
physeq_Cull <- subset_samples(physeq_prop, Timepoint == "Cull")

#----Plot recipient samples family bar graph----
# Define as family level - adjust for physeq object ***here*** as needed
barplot <- tax_glom(physeq_S0, taxrank = "Family")
barplot_recipient_family <- psmelt(barplot) #Melt phyloseq table into a data frame for plotting
#define unclassified
barplot_recipient_family$Family[is.na(barplot_donor_family$Family)] <- "Unclassified"

barplot <- ggplot(barplot_recipient_family, aes(x = Mouse_ID, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Group, scales = "free_x") +
  scale_fill_manual(values = family_colour) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(title = paste("Relative Abundance by group - S0"), x = "sample", y = "Relative Abundance")

barplot

#====Plotting at Genus level====
#Define taxonomy level
barplot_donor <- tax_glom(physeq_donor, taxrank = "Genus")
barplot_donor_genus <- psmelt(barplot_donor) #Melt physeq table into a data frame for plotting 

barplot_donor_genus$Genus[is.na(barplot_donor_genus$Genus)] <- "Unclassified"

#list of genus
unique_genus <-unique(taxa_filtered$Genus)
unique_genus

#generate a colour palette
genus_colour <- c(
  "Lactobacillus" = "red4",
  "Lactobacillaceae" = "red3",
  
  "Bacteroides" = "forestgreen",
  "Alistipes" = "darkolivegreen3",
  "Muribaculaceae" = "goldenrod2",
  
  "Akkermansia" = "dodgerblue3",
  
  "Prevotellaceae_UCG-001" = "blue4",
  
  "Lachnospiraceae_NK4A136_group" = "salmon2",
  "Lachnospiraceae_UCG-006" = "salmon3",
  "Lachnospiraceae_UCG-001" = "salmon4",
  "Lachnospiraceae" = "lightsalmon",
  
  "Roseburia" = "indianred2",
  "Blautia" = "rosybrown3",
  "Marvinbryantia" = "darksalmon",
  "[Eubacterium]_ventriosum_group" = "coral3",
  "[Eubacterium]_xylanophilum_group" = "coral4",
  "[Eubacterium]_coprostanoligenes_group" = "navy",
  "[Eubacterium]_brachy_group" = "tomato4",
  
  "Clostridium_sensu_stricto_1" = "tan4",
  "Clostridia_UCG-014" = "slateblue3",
  "Clostridia_vadinBB60_group" = "darkgoldenrod3",
  
  "Ruminococcus" = "#3C280D",
  "Oscillibacter" = "peachpuff3",
  "Oscillospiraceae" = "peachpuff1",
  
  "Faecalibaculum" = "goldenrod3",
  "Butyricicoccus" = "seagreen3",
  "Anaerotruncus" = "springgreen4",
  "Intestinimonas" = "mediumseagreen",
  
  "Turicibacter" = "orchid4",
  "Erysipelatoclostridium" = "deeppink3",
  "Erysipelotrichaceae" = "orchid3",
  "Erysipelatoclostridiaceae" = "deeppink4",
  
  "Bifidobacterium" = "aquamarine4",
  
  "Parasutterella" = "powderblue",
  
  "Escherichia-Shigella" = "lightskyblue4",
  "Enterococcus" = "violetred4",
  "Staphylococcus" = "firebrick3",
  "Staphylococcaceae" = "firebrick4",
  "Corynebacterium" = "darkseagreen",
  
  "Anaeroplasma" = "turquoise4",
  
  "Enterorhabdus" = "wheat4",
  "Coriobacteriaceae_UCG-002" = "burlywood4",
  
  "Tannerellaceae" = "palevioletred4",
  
  "Monoglobus" = "indianred3",
  "Colidextribacter" = "chocolate3",
  
  "Gastranaerophilales" = "olivedrab4",
  
  "Candidatus_Arthromitus" = "mediumpurple4",
  
  "ASF356" = "grey50",
  "HT002" = "grey60",
  "UCG-005" = "grey40",
  "UCG-010" = "mediumpurple4",
  "A2" = "grey70",
  "Incertae_Sedis" = "grey80",
  "Unclassified" = "grey90"
)

#----plot donor samples genus bar graph----

barplot_donor_genus <- ggplot(barplot_donor_genus, aes(x = Mouse_ID, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = genus_colour) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(title = paste("Relative Abundance (Genus level)"), x = "sample", y = "Relative Abundance")

barplot_donor_genus

#----plot recipients samples genus bar graph----
## Define as genus level - adjust for physeq object ***here*** as needed
barplot <- tax_glom(physeq_S0, taxrank = "Genus")
barplot_recipient_genus <- psmelt(barplot) #Melt phyloseq table into a data frame for plotting
#define unclassified
barplot_recipient_genus$Genus[is.na(barplot_recipient_genus$Genus)] <- "Unclassified"

##Plot genus bar graph
barplot_genus <- ggplot(barplot_recipient_genus, aes(x = Mouse_ID, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Group, scales = "free_x") +
  scale_fill_manual(values = genus_colour) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(title = paste("Relative Abundance by group - D15"), x = "sample", y = "Relative Abundance")

barplot_genus

#====Identify differential abundant genus in each donor (log fold change)====
#create a dataframe that contains **genus** from physeq_prop
genus_ps <- tax_glom(physeq_prop, taxrank = "Genus")
genus_df <- psmelt(genus_ps)
#create ordered time column
genus_df$Timepoint_ordered <- factor(
  genus_df$Timepoint,
  levels = c("S0","D0","D15","F7","F14","F21","Cull"),
  ordered = TRUE
)

#create a dataframe that contains **family** from physeq_prop
family_ps <- tax_glom(physeq_prop, taxrank = "Family")
family_df <-psmelt(family_ps)
#create ordered times
family_df$Timepoint_ordered <- factor(
  family_df$Timepoint,
  levels = c("S0","D0","D15","F7","F14","F21","Cull"),
  ordered = TRUE
)

#----Identify fold change of each genus by donor sample
#positive means enriched in donor-match FMT sample
#negative means enriched in auto-FMT sample
#look for abundance differences that are greater than e+02 for trajectory plots 

donor_genus_df <- genus_df %>%
  filter(Sample_ID %in% c("25-03302", "25-03303")) %>%
  select(Genus, Sample_ID, Abundance) %>%
  pivot_wider(names_from = Sample_ID, values_from = Abundance)

donor_genus_df <- donor_genus_df %>%
  mutate(log2FC = log2((`25-03303` + 1e-6) / (`25-03302` + 1e-6)))

#====longitudinal abundance of certain taxa over time====

#define taxa of interest ***Plug in taxa and run code from here***
taxa_of_interest <- "Muribaculaceae"

#CHOOSE:
#1. genus 
genus_target_df <- genus_df %>%
  filter(Genus == taxa_of_interest)

#OR

#2. family
family_target_df <- family_df %>%
  filter(Family == taxa_of_interest)

#1. define donor value 

donor_match_value <- family_target_df %>% #change to genus_target_df
  filter(Sample %in% c("25-03303")) %>%
  summarise(mean_abund = mean(Abundance)) %>%
  pull(mean_abund)

#1. plot genus trajectory

genus_trajetory <-ggplot(genus_target_df %>%
                           filter(Group %in% c("Donor match FMT", "AutoFMT", "Saline"), 
                                  Timepoint%in% c("S0","D0","D15","F7","F14","F21","Cull")),
       aes(x = Timepoint_ordered,
           y = Abundance,
           group = Mouse_ID,
           color = Group)) +
  geom_line(alpha = 0.3) +
  geom_point(alpha = 0.6) +
  stat_summary(aes(group = Group),
               fun = mean,
               geom = "line",
               linewidth = 1.5) +
  geom_hline(yintercept = donor_match_value, linetype = "dashed")+
  theme_classic() +
  labs(title = paste("Trajectory of", taxa_of_interest),
       y = "Relative abundance")

genus_trajetory

#2. plot family trajectory

family_trajetory <-ggplot(family_target_df %>%
                           filter(Group %in% c("Donor match FMT", "AutoFMT", "Saline"), 
                                  Timepoint%in% c("S0","D0","D15","F7","F14","F21","Cull")),
                         aes(x = Timepoint_ordered,
                             y = Abundance,
                             group = Mouse_ID,
                             color = Group)) +
  geom_line(alpha = 0.3) +
  geom_point(alpha = 0.6) +
  stat_summary(aes(group = Group),
               fun = mean,
               geom = "line",
               linewidth = 1.5) +
  geom_hline(yintercept = donor_match_value, linetype = "dashed")+
  theme_classic() +
  labs(title = paste("Family level trajectory of", taxa_of_interest),
       y = "Relative abundance")

family_trajetory

#====Saving graphs====
#Create consistent file path
save_plot <- function(plot_obj, filename, width = 12, height = 8) {
  ggsave(
    filename = file.path(results_dir, filename),
    plot = plot_obj,
    width = width,
    height = height
  )
}

#Saving graphs
save_plot(barplot_donor, "barplot_donor.pdf")
save_plot(barplot, "barplot_recipients.pdf")
save_plot(barplot_donor_genus, "barplot_donor_genus.pdf")
save_plot(barplot_genus, "barplot_genus.pdf", width = 16, height = 8)
save_plot(genus_trajetory, "Trajectory of Lachnospiraceae.pdf", width= 6, height=4)
save_plot(family_trajetory, "Family level trajectory of Muribaculaceae.pdf", width= 6, height=4)


