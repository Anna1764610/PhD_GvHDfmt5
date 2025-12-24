### Description: Analysing GvHDfmt5 16S dataset       ###
### Author: Anna Li, date: 17092025                   ###
### Email: anna.li@adelaide.edu.au                    ###

#Install phyloseq through Biocmanager command (provided by Bioconductor)
#BiocManager::install("phyloseq")

## Load Libraries
library(ape)              # Required for read.tree
library(phyloseq)         # Required for convert relative abundance to phyloseq object
library(vegan)            # Required for adonis2 - beta diversity
library(ggplot2)          # Required for graphing
library(ggpubr)           # Required for stat_compare_means when graphing shannon
library(gridExtra)        # Required for arranging plots (grid.arrange)
library(here)
library(readxl)

## Set up folders
project_dir <- here::here()

metadata_path <- file.path(
  project_dir,
  "Metadata",
  "GvHDfmt5_Metadata_R.xlsx"
)

fastq_dir <- 
  "~/OneDrive - University of Adelaide/PhD/Aim3_GvHDfmt study/GvHDfmt5/Raw_data/fastq/transfer_3397212_files_8dfb31e3"
#fastq is unnecessary to link to as QIIME2 has already crunched it but put as practice here for linking to cloud documents

##====Read in data====
raw_counts <- read.table(
  file.path(project_dir, "QIIME2_outputs", "feature-table.tsv"),
  sep = "\t",
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  comment.char = "",
  skip = 1
)

metadata <- readxl::read_excel(metadata_path)

tree = read.tree(
  file.path(project_dir,"QIIME2_outputs","tree.nwk"))
                 
taxa = read.delim(
  file.path(project_dir, "QIIME2_outputs", "taxonomy.tsv"),sep="\t",header=TRUE,row.names = 1)

#double check read in data appearance 
colnames(raw_counts) #Did R keep the sample ID?
head(taxa) #check on taxa
str(taxa) #what do the strings look like

#====Adjust taxa file====
# 1. Read the taxonomy file
taxa <- read.table(file.path(project_dir, "QIIME2_outputs","taxonomy.tsv"),
                   sep="\t",
                   header=TRUE,
                   stringsAsFactors=FALSE,
                   comment.char="")
# 2. Load needed libraries
library(tidyr)
library(dplyr)

# 3. Split Taxon column into ranks
taxa_clean <- taxa %>%
  separate(Taxon,
           into=c("Domain","Phylum","Class","Order","Family","Genus","Species"),
           sep=";", fill="right", remove=FALSE)
  
# 4. Remove prefixes (d__, p__, etc.)
taxa_clean <- taxa_clean %>%
  mutate(across(Domain:Species, ~ sub("^[a-z]__","", .)))

# 5. Drop the Taxon column
taxa_clean$Taxon <- NULL

# 6. set rownames to ASV ID (Feature.ID) rather than listing (1,2,3,4...)
rownames(taxa_clean) <- taxa_clean$Feature.ID

#check  
colnames(taxa_clean)
head(taxa_clean)

#====Assign empty, uncultured and others as unclassified====
# Replace na. and empty or uncultured strings everywhere with "Unclassified"
taxa_clean[is.na(taxa_clean)] <- "Unclassified"
taxa_clean[taxa_clean==""] <- "Unclassified"
taxa_clean[taxa_clean=="uncultured_bacterium"] <- "Unclassified"
taxa_clean[taxa_clean=="uncultured"] <- "Unclassified"
taxa_clean[taxa_clean=="Incertae_Sedis"] <- "Unclassified"
# Replace specific patterns in Species column
taxa_clean$Species[grepl("_sp\\.", taxa_clean$Species)] <- "Unclassified"
taxa_clean$Species[grepl("uncultured_", taxa_clean$Species)] <- "Unclassified"
taxa_clean$Species[grepl("mouse_gut", taxa_clean$Species)] <- "Unclassified"
taxa_clean$Species[grepl("unidentified", taxa_clean$Species)] <- "Unclassified"

# Optional: uncomment if you want to clean Genus column ##looks like they already are by QIIME2
# taxa_clean$Genus[grepl("GCA-", taxa_clean$Genus)] <- "Unclassified"

#====Cleaning data====
## Remove irrelevant data - SAGC negative samples
#remove rows where mouse ID is 'Negative'
metadata_clean<-metadata[metadata$Mouse_ID !="Negative",]
# Realign metadata and raw_counts
# First check names match...3 don't (since we removed Negatives)
colnames(raw_counts)
metadata_clean$Sample_ID
# Select only columns in raw_counts that are present in metadata_clean
raw_counts_clean <- raw_counts[, colnames(raw_counts) %in% metadata_clean$Sample_ID]
# now check numbers of samples align
colnames(raw_counts_clean)
metadata_clean$Sample_ID

#Add informative metadata by Mouse_ID 
# Load the Excel file
mouse_info <- read_excel(
  file.path(project_dir, "Metadata", "GvHDfmt5_study_parameter.xlsx"))
#Merge metadata files
metadata_merged <- merge(metadata_clean, mouse_info, by = "Mouse_ID", all.x = TRUE)
#remove column with same information from merging (group)
metadata_merged <- metadata_merged %>%
  select(-Group.y) %>%        # remove duplicate
  rename(Group = Group.x)     # rename cleanly

#remove more irrelevant data - mouse that died in procedure
#remove rows where HSCT_challange is 'EXCLUDE' (but not NA) #213 -> 209
metadata_clean2 <- metadata_merged[
  is.na(metadata_merged$HSCT_challenge) | metadata_merged$HSCT_challenge != "EXCLUDE", 
]

#cleanup file - restore rownames
rownames(metadata_clean2) <- metadata_clean2$Sample_ID

##Ensure 16S samples/raw_counts match new metadata file
#checks - to match, we need to end up with 209 samples as we have excluded based on metadata file
rownames(metadata_clean2) #209 rows of sample ID
colnames(raw_counts_clean) #213 col of sample ID

#Define overlapping samples
common_samples <- intersect(rownames(metadata_clean2), colnames(raw_counts_clean))

# Subset both objects to only the overlapping samples
raw_counts_clean2 <- raw_counts_clean[, common_samples]
metadata_clean2   <- metadata_clean2[common_samples, ]

colnames(raw_counts_clean2)

#====Tidying up the taxonomic data====
#To set anything that is unclassifeid at the genus level as it's family, and if still unclassified into order:
# Make sure Order and Genus are character columns
taxa_clean$Order <- as.character(taxa_clean$Order)
taxa_clean$Genus <- as.character(taxa_clean$Genus)
# Fill unclassified Genus with Family, then with Order
taxa_clean$Genus[taxa_clean$Genus == "Unclassified"] <- taxa_clean$Family[taxa_clean$Genus == "Unclassified"]
taxa_clean$Genus[taxa_clean$Genus == "Unclassified"] <- taxa_clean$Order[taxa_clean$Genus == "Unclassified"]

# Align raw_counts with taxa_clean (keep only features present in taxa)
# Make sure set rownames to ASV ID (Feature.ID) rather than listing (1,2,3,4...)
#rownames(taxa_clean) <- taxa_clean$Feature.ID
#raw_counts_clean2 <- raw_counts_clean2[rownames(taxa_clean), ]

#==== Filter out low-abundance taxa and low read samples ====

# Keep taxa with counts > 100 in at least 3 samples
raw_counts_filtered <- raw_counts_clean2[apply(raw_counts_clean2> 100, 1, sum) > 2, ]
  #1604 ASVs initially, 232 after filtering (object names: filtered) 

# Filter taxa_clean to match filtered raw_counts
taxa_filtered <- taxa_clean[rownames(raw_counts_filtered), ]

# Check read counts per sample in filtered samples - are there sufficient reads per sample?
sort(apply(raw_counts_filtered, 2, sum)) #sorts samples from least to greatest reads 
summary(apply(raw_counts_filtered, 2, sum))

# Total reads per sample
sample_depths <- colSums(raw_counts_filtered)

# Keep only samples with ≥5000 reads
keep_samples <- sample_depths >= 5000
raw_counts_keep <- raw_counts_filtered[, keep_samples]
metadata_keep   <- metadata_clean2[keep_samples, ]

  #209 samples initially, 207 after (object names: keep)

#====finding relative abundance (proportional counts) - can omit step ====

# Add total reads per sample in object metadata_keep
metadata_keep$read_num <- apply(raw_counts_keep, 2, sum)

# Convert counts to proportions (fraction)
counts_prop <- as.data.frame(prop.table(as.matrix(raw_counts_keep), 2))

#Extra filtering steps
# Keep taxa present in at least 1 sample (doesn't do much since already did in filtering step)
#counts_prop10 <- counts_prop[apply(counts_prop > 0, 1, sum) >= 1, ]
# Keep taxa with >0.1% relative abundance (0.001 in fraction) in more than 1 sample (no changes either)
#counts_prop10 <- counts_prop10[apply(counts_prop > 0.001, 1, sum) > 1, ]

##====Creating phyloseq object====

library(phyloseq)

# 1. Ensure filtered raw counts are properly formatted
# raw_counts_keep: ASVs x samples
OTU_raw <- otu_table(as.matrix(raw_counts_keep), taxa_are_rows = TRUE)

# 2. Ensure taxonomy is a matrix
tax_mat <- tax_table(as.matrix(taxa_filtered))

# 3. Make sure metadata rownames match sample IDs
rownames(metadata_keep) <- metadata_keep$Sample_ID
map_data <- sample_data(metadata_keep) #map_data is metadata_keep now recognised as phyloseq object 

# 4. Create phyloseq object with raw counts - use for alpha diversity, DESeq2 and other rarefraction-based metrics
physeq_raw <- phyloseq(OTU_raw, tax_mat, map_data)

# 5. Create a phyloseq object with relative abundances - use for bray-curtis beta diversity
physeq_prop <- transform_sample_counts(physeq_raw, function(x) x / sum(x))

# 6. Optional: check the object
physeq_raw
physeq_prop

# 7. Add tree to raw counts object
physeq_raw_tree <- merge_phyloseq(physeq_raw, phy_tree(tree))

# 8. Add tree to relative abundance object (optional)
# physeq_prop_tree <- merge_phyloseq(physeq_prop, phy_tree(tree))

#====Save RDS files: phyloseq objects====

saveRDS(physeq_raw, file=file.path(project_dir, "Phyloseq_objects", "physeq_raw.rds"))
saveRDS(physeq_prop, file=file.path(project_dir, "Phyloseq_objects", "physeq_prop.rds"))
saveRDS(metadata_keep, file=file.path(project_dir, "Phyloseq_objects", "metadata_keep.rds"))
saveRDS(physeq_raw_tree, file=file.path(project_dir, "Phyloseq_objects", "physeq_raw_tree.rds"))

