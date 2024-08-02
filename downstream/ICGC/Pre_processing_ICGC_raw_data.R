
library(tidyverse)
library(ggplot2)

location <- "local" # 'local' or 'hpc'

path_SSMs <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/genomics/raw_ICGC/simple_somatic_mutation.open.tsv")
path_SSMs_with_AF <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/genomics/raw_ICGC/simple_somatic_mutation.open.with_AF.tsv")
path_STSMs <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/genomics/raw_ICGC/structural_somatic_mutation.tsv")

path_results <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/genomics/pre_processed_ICGC/")


if(location == "hpc"){
  path_SSMs <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/data/genomics/raw_ICGC/simple_somatic_mutation.open.tsv")
  path_SSMs_with_AF <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/data/genomics/raw_ICGC/simple_somatic_mutation.open.with_AF.tsv")
  path_STSMs <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/data/genomics/raw_ICGC/structural_somatic_mutation.tsv")
  
  path_results <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/data/genomics/pre_processed_ICGC/")
}
  

### SIMPLE SOMATIC MUTATIONS ###

# Raw ICGC file contains 17'889'165 instances
SSMs <- read_tsv(path_SSMs)
dim(SSMs)


# Multiple variant callers were used to call variants
table(SSMs[,15])

# But most calls are equal among callers. SSMs_distinct contains 5'012'977 fields
SSMs_distinct <- distinct(SSMs[, -15])
dim(SSMs_distinct)


## Among all the SSMs, some are present only once, others multiple times

# Mutations in which, for each icgc_mutation_id, only one instance is present
SSMs_matching <- SSMs_distinct %>% group_by(icgc_mutation_id) %>%
  mutate(., n_ist = length(icgc_mutation_id)) %>%
  dplyr::filter(., n_ist == 1) %>%
  ungroup()

# Mutations for which, for each icgc_mutation_id, more istances are present
SSMs_not_matching <- SSMs_distinct %>% group_by(icgc_mutation_id) %>%
  mutate(., n_ist = length(icgc_mutation_id)) %>%
  dplyr::filter(., n_ist != 1) %>%
  ungroup()


# Multiple instances of the same icgc_mutation_id differ for one or few of the other fields (such as: icgc_donor_id, icgc_sample_id, etc)
# I.e. because the SAME mutation (with a unique icgc_mutation_id) is found in MULTIPLE donors (several icgc_donor_id)
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


## Re-joining with matching SSMs
SSMs_new <- rbind(SSMs_matching, SSMs_not_matching)

#write_tsv(SSMs_new, fs::path(path_results, "simple_somatic_mutation.open.matching_calls.tsv"))




### STRUCTURAL SOMATIC MUTATIONS ###

# STSMs 
STSMs <- read_tsv(path_STSMs) # it contains by default only unique calls
STSMs <- STSMs[, c(1:23)]

# Each STSM has a 'from' breakpoint and a 'to' breakpoint. Split the two info to maintain both as separate bkpts
STSMs_coords <- STSMs
STSMs_coords_from <- STSMs_coords[, c(1:16)]
STSMs_coords_from$class <- rep("from", dim(STSMs_coords_from)[1])
colnames(STSMs_coords_from) <- c(colnames(STSMs_coords_from)[1:11], c("chrom", "chrom_bkpt", "strand", "range", "chrom_flanking_seq", "class"))

STSMs_coords_to <- STSMs_coords[, c(1:11, 17:21)]
STSMs_coords_to$class <- rep("to", dim(STSMs_coords_to)[1])
colnames(STSMs_coords_to) <- c(colnames(STSMs_coords_to)[1:11], c("chrom", "chrom_bkpt", "strand", "range", "chrom_flanking_seq", "class"))

STSMs_coords_whole <- rbind(
  STSMs_coords_from, STSMs_coords_to
)

STSMs_coords_whole$sv_id_x_class <- 
  paste(STSMs_coords_whole$sv_id, STSMs_coords_whole$class, sep = "_")
STSMs_coords_whole <- relocate(STSMs_coords_whole, sv_id_x_class, .after = sv_id)

#write_tsv(STSMs_coords_whole, fs::path(path_results, "structural_somatic_mutation.preprocessed.tsv"))



### SIMPLE SOMATIC MUTATIONS WITH AF ###

# Raw ICGC file contains 17'889'165 instances
SSMs <- read_tsv(path_SSMs_with_AF)
dim(SSMs)


# Multiple variant callers were used to call variants
table(SSMs[,16])

# But most calls are equal among callers. SSMs_distinct contains 5'012'977 fields
SSMs_distinct <- distinct(SSMs[, -16])
dim(SSMs_distinct)


## Among all the SSMs, some are present only once, others multiple times

# Mutations in which, for each icgc_mutation_id, only one instance is present
SSMs_matching <- SSMs_distinct %>% group_by(icgc_mutation_id) %>%
  mutate(., n_ist = length(icgc_mutation_id)) %>%
  dplyr::filter(., n_ist == 1) %>%
  ungroup()

# Mutations for which, for each icgc_mutation_id, more istances are present
SSMs_not_matching <- SSMs_distinct %>% group_by(icgc_mutation_id) %>%
  mutate(., n_ist = length(icgc_mutation_id)) %>%
  dplyr::filter(., n_ist != 1) %>%
  ungroup()


# Multiple instances of the same icgc_mutation_id differ for one or few of the other fields (such as: icgc_donor_id, icgc_sample_id, etc)
# I.e. because the SAME mutation (with a unique icgc_mutation_id) is found in MULTIPLE donors (several icgc_donor_id)
SSMs_var_sum <- SSMs_not_matching %>% group_by(icgc_mutation_id) %>%
  dplyr::summarise(
    n_donors = length(unique(icgc_donor_id)), 
    n_specimen = length(unique(icgc_specimen_id)), 
    n_sample = length(unique(icgc_sample_id)), 
    n_project = length(unique(project_code)), 
    n_VAFs = length(unique(total_read_count)))
summary(SSMs_var_sum)


# ----- #
# Filter out variants WITHOUT total_read_count - since we need AFs, we are not interested in those not having this info
SSMs_not_matching_with_AF <- SSMs_not_matching %>% filter(., !is.na(total_read_count))
# Some mutations are called in multiple samples (for the same donor) and have different total_read_count / mutant_allele_read_count
# We select the mutation having the highest total_read_count
SSMs_not_matching_with_AF_clean <- SSMs_not_matching_with_AF %>% group_by(icgc_mutation_id, icgc_donor_id) %>% 
  dplyr::filter(total_read_count == max(total_read_count)) %>% 
  ungroup()
# Some variants have equal total_read_counts - slect the first copy only (random)
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

## Re-joining with matching SSMs
SSMs_new <- rbind(SSMs_matching, SSMs_not_matching)
SSMs_new$AF <- round(SSMs_new$mutant_allele_read_count / SSMs_new$total_read_count, 3)
  
#write_tsv(SSMs_new, fs::path(path_results, "simple_somatic_mutation.open.matching_calls.with_AFs.tsv"))


