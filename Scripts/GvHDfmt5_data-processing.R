### Description: Analysing GvHDfmt5 16S dataset       ###
### Author: Anna Li, date: 17092025                   ###
### Email: anna.li@adelaide.edu.au                    ###


## Load Libraries
library(ape)              # Required for read.tree
library(phyloseq)         # Required for convert relative abundance to phyloseq object
library(vegan)            # Required for adonis2 - beta diversity
library(ggplot2)          # Required for graphing
library(ggpubr)           # Required for stat_compare_means when graphing shannon
library(gridExtra)        # Required for arranging plots (grid.arrange)
library(here)
library(readxl)
library(tidyr)
library(dplyr)

## Set up folders
project_dir <- here::here()

metadata_path <- file.path(
  project_dir,
  "Metadata",
  "GvHDfmt5_Metadata_R.xlsx"
)

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

#----Adjust taxa file----
#Split Taxon column into ranks
taxa_clean <- taxa %>%
  separate(Taxon,
           into=c("Domain","Phylum","Class","Order","Family","Genus","Species"),
           sep=";", fill="right", remove=FALSE)%>% mutate(across(Domain:Species, ~ sub("^[a-z]__","", .)))

#Drop the Taxon column (after checking, otherwise can get rid of "remove = FALSE" in prior code)
taxa_clean$Taxon <- NULL

# Replace na. and empty or uncultured strings everywhere with "Unclassified"
taxa_clean[is.na(taxa_clean)] <- "Unclassified"
taxa_clean[taxa_clean==""] <- "Unclassified"
taxa_clean[taxa_clean=="uncultured_bacterium"] <- "Unclassified"
taxa_clean[taxa_clean=="uncultured"] <- "Unclassified"

#To set anything that is unclassified at the genus level as its family, and if still unclassified into order:
# Make sure Order and Genus are character columns
taxa_clean$Order <- as.character(taxa_clean$Order)
taxa_clean$Genus <- as.character(taxa_clean$Genus)
# Fill unclassified Genus with Family, then with Order
taxa_clean$Genus[taxa_clean$Genus == "Unclassified"] <- taxa_clean$Family[taxa_clean$Genus == "Unclassified"]
taxa_clean$Genus[taxa_clean$Genus == "Unclassified"] <- taxa_clean$Order[taxa_clean$Genus == "Unclassified"]

#----Clean data----
## Remove irrelevant data - SAGC negative samples
#remove SAGC negative from metadata
metadata_clean<-metadata[metadata$Mouse_ID !="Negative",]
## Realign metadata and raw_counts
# First check names match...3 don't (since we removed Negatives)
colnames(raw_counts)
metadata_clean$Sample_ID
# Select only columns in raw_counts that are present in metadata_clean
raw_counts_clean <- raw_counts[, colnames(raw_counts) %in% metadata_clean$Sample_ID]
# now check numbers of samples align
colnames(raw_counts_clean)
metadata_clean$Sample_ID

##Add informative metadata by Mouse_ID
# Load the Excel file
mouse_info <- read_excel(
  file.path(project_dir, "Metadata", "GvHDfmt5_study_parameter.xlsx"))
#Merge metadata files
metadata_merged <- merge(metadata_clean, mouse_info, by = "Mouse_ID", all.x = TRUE)
#remove column with same information from merging (group)
metadata_merged <- metadata_merged %>%
  select(-Group.y) %>%        # remove duplicate
  rename(Group = Group.x)     # rename cleanly

##remove more irrelevant data - mouse that died in procedure
#remove rows where HSCT_challenge is 'EXCLUDE' (but not NA) #213 obs -> 209 obs
metadata_clean2 <- metadata_merged[
  is.na(metadata_merged$HSCT_challenge) | metadata_merged$HSCT_challenge != "EXCLUDE", 
]
#cleanup file - restore SampleID as rownames 
rownames(metadata_clean2) <- metadata_clean2$Sample_ID

##Ensure 16S samples/raw_counts match new metadata file
#checks - to match, we need to end up with 209 samples as we have excluded based on metadata file
rownames(metadata_clean2) #209 rows of sample ID
colnames(raw_counts_clean) #213 col of sample ID
#need to remove 4 cols from raw_counts_clean

#Define overlapping samples
common_samples <- intersect(rownames(metadata_clean2), colnames(raw_counts_clean))

# Subset both objects to only the overlapping samples
raw_counts_clean2 <- raw_counts_clean[, common_samples]
metadata_clean2   <- metadata_clean2[common_samples, ]

colnames(raw_counts_clean2)

#==== Filter out low-abundance taxa and low read samples ====
# Keep taxa with counts > 100 in at least 2 samples (TWO SAMPLES NOT 2 - made adjustment)
raw_counts_filtered <- raw_counts_clean2[apply(raw_counts_clean> 100, 1, sum) > 1, ]
  #1604 ASVs initially, 257 after filtering (object names: filtered) 

# Filter taxa_clean to match filtered raw_counts
taxa_filtered <- taxa_clean[rownames(raw_counts_filtered), ]

# Check read counts per sample in filtered samples - are there sufficient reads per sample?
sort(apply(raw_counts_filtered, 2, sum)) #sorts samples from least to greatest reads 
summary(apply(raw_counts_filtered, 2, sum)) #gives summary stats

# Keep only samples with ≥5000 reads - should remove two samples from the analysis where read depth is too low
sample_depths <- colSums(raw_counts_filtered)
keep_samples <- sample_depths >= 5000
raw_counts_keep <- raw_counts_filtered[, keep_samples]
metadata_keep   <- metadata_clean2[keep_samples, ] #realign metadata

  #209 samples initially, 207 after (object names: keep)

##====Creating phyloseq object====

## Create phyloseq object components recognised by phyloseq command
# Create OTU_raw from raw_counts_keep as a matrix
OTU_raw <- otu_table(as.matrix(raw_counts_keep), taxa_are_rows = TRUE)

# Create tax_mat file from taxa_filtered as a matrix
tax_mat <- tax_table(as.matrix(taxa_filtered))

# Create map_data
# Make sure metadata rownames match sample IDs
rownames(metadata_keep) <- metadata_keep$Sample_ID
map_data <- sample_data(metadata_keep) #map_data is metadata_keep now recognised as phyloseq object

## Create phyloseq object with raw counts - use for alpha diversity, DESeq2 and other rarefraction-based metrics
physeq_raw <- phyloseq(OTU_raw, tax_mat, map_data)

# Create a phyloseq object with relative abundances - use for bray-curtis beta diversity
physeq_prop <- transform_sample_counts(physeq_raw, function(x) x / sum(x))

# Optional: check the object has converted into relative abundance
sample_sums(physeq_raw)[1:5]
sample_sums(physeq_prop)[1:5] #should all =1

# Add tree to raw counts object
physeq_raw_tree <- merge_phyloseq(physeq_raw, phy_tree(tree))

# Add tree to relative abundance object (optional)
# physeq_prop_tree <- merge_phyloseq(physeq_prop, phy_tree(tree))

#====Save RDS files: phyloseq objects====
#saveRDS(physeq_raw, file=file.path(project_dir, "Phyloseq_objects", "physeq_raw.rds"))
#saveRDS(physeq_prop, file=file.path(project_dir, "Phyloseq_objects", "physeq_prop.rds"))
#saveRDS(map_data, file=file.path(project_dir, "Phyloseq_objects", "map_data.rds"))
#saveRDS(physeq_raw_tree, file=file.path(project_dir, "Phyloseq_objects", "physeq_raw_tree.rds"))


