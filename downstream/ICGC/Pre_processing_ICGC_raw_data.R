
library(tidyverse)
library(ggplot2)

location <- "local" # 'local' or 'hpc'

# Define file paths for local environment
path_SSMs <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/genomics/raw_ICGC/simple_somatic_mutation.open.tsv")
path_SSMs_with_AF <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/genomics/raw_ICGC/simple_somatic_mutation.open.with_AF.tsv")
path_STSMs <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/genomics/raw_ICGC/structural_somatic_mutation.tsv")

path_results <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/genomics/pre_processed_ICGC/")

# Override file paths if using HPC environment
if(location == "hpc"){
  path_SSMs <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/data/genomics/raw_ICGC/simple_somatic_mutation.open.tsv")
  path_SSMs_with_AF <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/data/genomics/raw_ICGC/simple_somatic_mutation.open.with_AF.tsv")
  path_STSMs <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/data/genomics/raw_ICGC/structural_somatic_mutation.tsv")
  
  path_results <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/data/genomics/pre_processed_ICGC/")
}
  

##


### SIMPLE SOMATIC MUTATIONS ###


# Load the raw ICGC simple somatic mutation (SSM) data
# The dataset contains 17,889,165 entries
SSMs <- read_tsv(path_SSMs)
dim(SSMs)


# Check how many distinct variant callers (column 15) are used in the data
table(SSMs[,15])

# Remove duplicate entries based on all columns except the variant caller
# Most calls are equal among callers. The resulting dataset contains 5,012,977 entries
SSMs_distinct <- distinct(SSMs[, -15])
dim(SSMs_distinct)


## Among all the SSMs, some are present only once, others multiple times

# Select mutations that are unique, i.e., they appear only once per icgc_mutation_id
SSMs_matching <- SSMs_distinct %>% group_by(icgc_mutation_id) %>%
  mutate(., n_ist = length(icgc_mutation_id)) %>%
  dplyr::filter(., n_ist == 1) %>%
  ungroup()

# Select mutations that have multiple instances per icgc_mutation_id
SSMs_not_matching <- SSMs_distinct %>% group_by(icgc_mutation_id) %>%
  mutate(., n_ist = length(icgc_mutation_id)) %>%
  dplyr::filter(., n_ist != 1) %>%
  ungroup()


# Multiple instances of the same icgc_mutation_id differ for one or few of the other fields (such as: icgc_donor_id, icgc_sample_id, etc)
# This occurs because the SAME mutation (with a unique icgc_mutation_id) is found in MULTIPLE donors (several icgc_donor_id)
SSMs_var_sum <- SSMs_not_matching %>% group_by(icgc_mutation_id) %>%
  dplyr::summarise(
    n_donors = length(unique(icgc_donor_id)), 
    n_specimen = length(unique(icgc_specimen_id)), 
    n_sample = length(unique(icgc_sample_id)), 
    n_project = length(unique(project_code)))
summary(SSMs_var_sum)


# For now, we are not interested in keeping mutations that are present in multiple instances because they differ for the sample_id
# or specimen_id, BUT come from the same donor
SSMs_not_matching$n_ist <- NULL
SSMs_not_matching$distinct_mut <- !duplicated(SSMs_not_matching[, -c(3:5)]) # removing ONLY duplicated mutations equal in ALL fields except for project_code, sample, specimen
sum(!SSMs_not_matching$distinct_mut) # These are 4'068 duplicates

# Remove them
SSMs_not_matching <- SSMs_not_matching[SSMs_not_matching$distinct_mut == T, ]
SSMs_not_matching$n_ist <- NULL
SSMs_not_matching$distinct_mut <- NULL

SSMs_matching$n_ist <- NULL


# Combine the unique and non-duplicated SSMs into a single dataset
SSMs_new <- rbind(SSMs_matching, SSMs_not_matching)

#write_tsv(SSMs_new, fs::path(path_results, "simple_somatic_mutation.open.matching_calls.tsv"))


##


### SIMPLE SOMATIC MUTATIONS WITH AF ###

# Load the raw ICGC simple somatic mutation (SSM) data with allele frequency (AF) information
# The dataset contains 17,889,165 entries
SSMs <- read_tsv(path_SSMs_with_AF)
dim(SSMs)


# Check how many distinct variant callers (column 16) are used in the data
table(SSMs[,16])

# Remove duplicate entries based on all columns except the variant caller
# The resulting dataset contains 5,012,977 entries
SSMs_distinct <- distinct(SSMs[, -16])
dim(SSMs_distinct)


## Among all the SSMs, some are present only once, others multiple times

# Select mutations that are unique, i.e., they appear only once per icgc_mutation_id
SSMs_matching <- SSMs_distinct %>% group_by(icgc_mutation_id) %>%
  mutate(., n_ist = length(icgc_mutation_id)) %>%
  dplyr::filter(., n_ist == 1) %>%
  ungroup()

# Select mutations that have multiple instances per icgc_mutation_id
SSMs_not_matching <- SSMs_distinct %>% group_by(icgc_mutation_id) %>%
  mutate(., n_ist = length(icgc_mutation_id)) %>%
  dplyr::filter(., n_ist != 1) %>%
  ungroup()


# Multiple instances of the same icgc_mutation_id differ for one or few of the other fields (such as: icgc_donor_id, icgc_sample_id, etc)
# This occurs because the SAME mutation (with a unique icgc_mutation_id) is found in MULTIPLE donors (several icgc_donor_id)
SSMs_var_sum <- SSMs_not_matching %>% group_by(icgc_mutation_id) %>%
  dplyr::summarise(
    n_donors = length(unique(icgc_donor_id)), 
    n_specimen = length(unique(icgc_specimen_id)), 
    n_sample = length(unique(icgc_sample_id)), 
    n_project = length(unique(project_code)), 
    n_VAFs = length(unique(total_read_count)))
summary(SSMs_var_sum)


# ----- #

# Filter out variants that do not have total_read_count information
# Since we need AFs, we are not interested in mutations lacking this info
SSMs_not_matching_with_AF <- SSMs_not_matching %>% filter(., !is.na(total_read_count))

# For mutations found in multiple samples for the same donor, select the mutation with the highest total_read_count
SSMs_not_matching_with_AF_clean <- SSMs_not_matching_with_AF %>% group_by(icgc_mutation_id, icgc_donor_id) %>% 
  dplyr::filter(total_read_count == max(total_read_count)) %>% 
  ungroup()

# For variants with equal total_read_count, select the first copy (random selection)
SSMs_not_matching_with_AF %>% group_by(icgc_mutation_id, icgc_donor_id) %>% 
  slice(1)

# ---- #


# For now, we are not interested in keeping mutations that are present in multiple instances because they differ for the sample_id
# or specimen_id, BUT come from the same donor
# THIS WILL KEEP DUPLICATED VARIANTS WITH DIFFERENT AFs - WILL CONSIDER ICGC_SAMPLE_ID INSTEAD OF ICGC_DONOR_ID
SSMs_not_matching$n_ist <- NULL
SSMs_not_matching$distinct_mut <- !duplicated(SSMs_not_matching[, -c(3:5)]) # removing ONLY duplicated mutations equal in ALL fields except for project_code, sample, specimen
sum(!SSMs_not_matching$distinct_mut) # These are 4'068 duplicates

# Remove them
SSMs_not_matching <- SSMs_not_matching[SSMs_not_matching$distinct_mut == T, ]
SSMs_not_matching$n_ist <- NULL
SSMs_not_matching$distinct_mut <- NULL

SSMs_matching$n_ist <- NULL

# Combine the unique and non-duplicated SSMs into a single dataset
SSMs_new <- rbind(SSMs_matching, SSMs_not_matching)

# Calculate the allele frequency (AF) for each mutation and add it as a new column
SSMs_new$AF <- round(SSMs_new$mutant_allele_read_count / SSMs_new$total_read_count, 3)
  
# Save the cleaned SSM dataset with AFs
#write_tsv(SSMs_new, fs::path(path_results, "simple_somatic_mutation.open.matching_calls.with_AFs.tsv"))


##


### STRUCTURAL SOMATIC MUTATIONS ###


# Load the raw ICGC structural somatic mutation (STSM) data
# The dataset contains only unique calls by default
STSMs <- read_tsv(path_STSMs) 
STSMs <- STSMs[, c(1:23)]

# Each STSM has a 'from' breakpoint and a 'to' breakpoint
# Split the dataset into two: one for 'from' breakpoints and one for 'to' breakpoints
STSMs_coords <- STSMs
STSMs_coords_from <- STSMs_coords[, c(1:16)]
STSMs_coords_from$class <- rep("from", dim(STSMs_coords_from)[1])
colnames(STSMs_coords_from) <- c(colnames(STSMs_coords_from)[1:11], c("chrom", "chrom_bkpt", "strand", "range", "chrom_flanking_seq", "class"))

STSMs_coords_to <- STSMs_coords[, c(1:11, 17:21)]
STSMs_coords_to$class <- rep("to", dim(STSMs_coords_to)[1])
colnames(STSMs_coords_to) <- c(colnames(STSMs_coords_to)[1:11], c("chrom", "chrom_bkpt", "strand", "range", "chrom_flanking_seq", "class"))

# Combine the 'from' and 'to' breakpoints into a single dataset
STSMs_coords_whole <- rbind(
  STSMs_coords_from, STSMs_coords_to
)

# Create a new column that combines the sv_id and class for each entry
STSMs_coords_whole$sv_id_x_class <- 
  paste(STSMs_coords_whole$sv_id, STSMs_coords_whole$class, sep = "_")
STSMs_coords_whole <- relocate(STSMs_coords_whole, sv_id_x_class, .after = sv_id)

# Save the preprocessed STSM dataset
#write_tsv(STSMs_coords_whole, fs::path(path_results, "structural_somatic_mutation.preprocessed.tsv"))


##