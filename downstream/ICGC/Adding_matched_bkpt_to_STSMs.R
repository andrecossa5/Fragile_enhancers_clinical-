
library(tidyverse)  # Load the tidyverse library for data manipulation

# Read in the data for structural somatic mutations (STSMs) overlapping with CTIP and GRHL enhancers
STSMs_ctip <- read_tsv("/Users/ieo6983/Desktop/enhancers_project_old/per_ivan/STSMs_overlapping_CTIP_enhancers.tsv")
STSMs_grhl <- read_tsv("/Users/ieo6983/Desktop/enhancers_project_old/per_ivan/STSMs_overlapping_GRHL_enhancers.tsv")

# Read in the preprocessed and raw ICGC STSMs data
STSMs_all <- read_tsv("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/genomics/pre_processed_ICGC/structural_somatic_mutation.preprocessed.tsv")
STSMs_raw <- read_tsv("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/genomics/raw_ICGC/structural_somatic_mutation.tsv")


##


# Prepare STSMs data for merging
STSMs_all_to_merge <- STSMs_all %>% 
  dplyr::select(., chrom, chrom_bkpt, icgc_donor_id, variant_type, sv_id, class)  # Select relevant columns
colnames(STSMs_all_to_merge)[1:2] <- c("seqnames_sbj2", "start_sbj2")  # Rename columns for consistency
STSMs_all_to_merge <- STSMs_all_to_merge %>% 
  mutate(., end_sbj2 = .$start_sbj2, .after = start_sbj2)  # Add an "end_sbj2" column identical to "start_sbj2"
STSMs_all_to_merge$seqnames_sbj2 <- paste0("chr", STSMs_all_to_merge$seqnames_sbj2)  # Prefix chromosome names with "chr"

# Check if all sv_ids in STSMs_ctip and STSMs_grhl are present in STSMs_all_to_merge
sum(!STSMs_ctip$sv_id %in% STSMs_all_to_merge$sv_id)  # Count sv_ids in STSMs_ctip not found in STSMs_all_to_merge
sum(!STSMs_grhl$sv_id %in% STSMs_all_to_merge$sv_id)  # Count sv_ids in STSMs_grhl not found in STSMs_all_to_merge


##


## Add second breakpoint information (bkpt) in two steps

# Step 1: Get unique breakpoints from the overlaps dataframe
STSMs_ctip1 <- STSMs_ctip[!duplicated(STSMs_ctip$sv_id), ]  # Select unique sv_ids from STSMs_ctip
STSMs_ctip2 <- STSMs_ctip[duplicated(STSMs_ctip$sv_id), ]  # Select duplicated sv_ids from STSMs_ctip

# Step 2: Prepare data for merging with the second breakpoint information
ctip_to_merge <- STSMs_all_to_merge %>% 
  semi_join(., STSMs_ctip, by = "sv_id") %>%  # Keep only rows from STSMs_all_to_merge with matching sv_ids in STSMs_ctip
  anti_join(., STSMs_ctip, by = c("sv_id", "class"))  # Remove rows with matching sv_id and class in STSMs_ctip

left_join()  # Placeholder for potential merging (appears incomplete in the original script)


##


# Process STSMs data for CTIP enhancers

# Initialize an empty dataframe to store results
STSMs_ctip_new <- data.frame()

# Iterate over each row in STSMs_ctip
for(i in 1:dim(STSMs_ctip)[1]){
  row_1 <- STSMs_ctip[i, ]  # Extract the current row from STSMs_ctip
  row_2 <- STSMs_all_to_merge[ (STSMs_all_to_merge$sv_id == row_1$sv_id & !STSMs_all_to_merge$class == row_1$class) , ]  # Find the corresponding row in STSMs_all_to_merge with the same sv_id but different class
  connected_rows <- cbind(row_1, row_2)  # Combine the two rows side by side
  
  STSMs_ctip_new <- rbind(STSMs_ctip_new, connected_rows)  # Append the combined row to the new dataframe
}
colnames(STSMs_ctip_new) <- c(
  colnames(STSMs_ctip_new)[1:16], 
  c("icgc_donor_id2", "variant_type2", "sv_id2", "class2") # Rename the columns for the second breakpoint information
)

# Check for consistency between sv_id and sv_id2, and between class and class2
sum(!STSMs_ctip_new$sv_id == STSMs_ctip_new$sv_id2)  # Count rows where sv_id does not match sv_id2
sum(STSMs_ctip_new$class == STSMs_ctip_new$class2)  # Count rows where class matches class2


# Process STSMs data for GRHL2 enhancers

# Initialize an empty dataframe to store results
STSMs_grhl_new <- data.frame()

# Iterate over each row in STSMs_grhl
for(i in 1:dim(STSMs_grhl)[1]){
  row_1 <- STSMs_grhl[i, ]  # Extract the current row from STSMs_grhl
  row_2 <- STSMs_all_to_merge[ (STSMs_all_to_merge$sv_id == row_1$sv_id & !STSMs_all_to_merge$class == row_1$class) , ]  # Find the corresponding row in STSMs_all_to_merge with the same sv_id but different class
  connected_rows <- cbind(row_1, row_2)  # Combine the two rows side by side
  
  STSMs_grhl_new <- rbind(STSMs_grhl_new, connected_rows)  # Append the combined row to the new dataframe
}
colnames(STSMs_grhl_new) <- c(
  colnames(STSMs_grhl_new)[1:16], 
  c("icgc_donor_id2", "variant_type2", "sv_id2", "class2") # Rename the columns for the second breakpoint information
)

# Check for consistency between sv_id and sv_id2, and between class and class2
sum(!STSMs_grhl_new$sv_id == STSMs_grhl_new$sv_id2)  # Count rows where sv_id does not match sv_id2
sum(STSMs_grhl_new$class == STSMs_grhl_new$class2)  # Count rows where class matches class2


##


# The following lines were used to write the results to files
#write_tsv(STSMs_ctip_new, "/Users/ieo6983/Desktop/enhancers_project_old/per_ivan/STSMs_overlapping_CTIP_enhancers.both_bkpts.tsv")
#write_tsv(STSMs_grhl_new, "/Users/ieo6983/Desktop/enhancers_project_old/per_ivan/STSMs_overlapping_GRHL_enhancers.both_bkpts.tsv")
