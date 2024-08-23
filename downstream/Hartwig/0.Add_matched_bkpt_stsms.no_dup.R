
library(tidyverse)  # Load the tidyverse package for data manipulation

# Set a seed for reproducibility
SEED <- 4321
set.seed(SEED)

# Flag to determine if the output should be saved
save_output <- T

# Define file paths for input and output data
path_stsms_all <- fs::path("/hpcnfs/scratch/P_PGP_FRAGILE_ENHANCERS/results/data/Hartwig_all_stsms_info.tsv")
path_results <- fs::path("/hpcnfs/scratch/P_PGP_FRAGILE_ENHANCERS/results/NEW/data/")
  

##


# Read in the structural variant data from the specified file path
stsms_all <- read_tsv(path_stsms_all)

# Each row in the dataset should contain only 1 breakpoint
print("Each row contains only 1 bkpt")
print(sum(!stsms_all$start == stsms_all$end) == 0)  # Check if start and end positions are the same for all rows

# Display the total number of entries (single breakpoints)
print("Total number of entries - single bkpts: ")
print(dim(stsms_all)[1])  # 871'123 total rows

# Display the total number of unique IDs and PARIDs
print("Total number of unique IDs and PARIDs:")
print(length(unique(stsms_all$ID)))  # 836'226 unique IDs
print(length(unique(stsms_all$PARID)))  # 789'061 unique PARIDs

# Check for NA values in IDs and PARIDs
print("Total number of unique IDs and PARIDs:")
sum(is.na(stsms_all$ID))  # Count of NAs in ID column
sum(is.na(stsms_all$PARID))  # Count of NAs in PARID column

# Check for duplicated IDs and PARIDs
print("Duplicates in IDs and PARIDs:")
print(sum(duplicated(stsms_all$ID)))  # Count of duplicated IDs
print(sum(duplicated(stsms_all$PARID)))  # Count of duplicated PARIDs

#' @Notes
#' Duplicated IDs are present, and can be associated to the SAME PARID, even though the variant type is different
#' example:
#' TYPE: INS  ID: gridss0_12556h  PARID: gridss0_12556o
#' TYPE: INV  ID: gridss0_12556h  PARID: gridss0_12556o
#' 
#' TYPE: INS  ID: gridss0_12556o  PARID: gridss0_12556h
#' TYPE: INV  ID: gridss0_12556o  PARID: gridss0_12556h
#' 
#' I don't know why this ambiguity is present, 
#' could be due to duplicated ID - PARID entries, called for different SAMPLEs
#' but it looks like the only solution is to check both ID - PARID match and TYPE - TYPE match


##


# Initialize an empty data frame to store the final results
STSMs_all_new <- data.frame()
# Initialize a vector to keep track of rows that should be skipped to avoid duplication
i_to_skip <- c()

# Loop over each row in the input data
for(i in 1:dim(stsms_all)[1]){
  print(i)  # Print the current iteration number
  
  # Skip rows that have already been processed
  if(i %in% i_to_skip){
    next
  }
  
  # Skip rows that are entirely empty
  if( sum(!is.na(stsms_all[i, ])) == 0 ){
    next
  }
  
  # Extract the current row, which contains ID and PARID information
  row_1 <- stsms_all[i, ] 
  
  # If PARID is missing, create a placeholder row with NA values
  if(is.na(row_1$PARID)){
    row_2 <- colnames(stsms_all)
    df <- data.frame(matrix(NA, nrow=1, ncol=length(row_2)))
    colnames(df) <- row_2
    connected_rows <- cbind(row_1, df)  # Combine the original row with the placeholder
    
    # Append the combined row to the new dataset
    STSMs_all_new <- rbind(STSMs_all_new, connected_rows)
  } 
  else {
    # Find the row containing the PARID breakpoint, matching by TYPE and SAMPLE
    row_2 <- stsms_all[(stsms_all$ID == row_1$PARID & stsms_all$TYPE == row_1$TYPE & stsms_all$SAMPLE == row_1$SAMPLE), ]
    row_2_i <- which(stsms_all$ID == row_1$PARID & stsms_all$TYPE == row_1$TYPE & stsms_all$SAMPLE == row_1$SAMPLE)  # Get the index of the matching row
    
    # If the PARID row is not present (possibly filtered out), create a placeholder row with NA values
    if(dim(row_2)[1] == 0){
      row_2 <- colnames(stsms_all)
      df <- data.frame(matrix(NA, nrow=1, ncol=length(row_2)))
      colnames(df) <- row_2
      connected_rows <- cbind(row_1, df)  # Combine the original row with the placeholder
      
      # Append the combined row to the new dataset
      STSMs_all_new <- rbind(STSMs_all_new, connected_rows)
    }
    # If the PARID row is present, combine it with the original row
    else {
      # Mark the PARID row index to be skipped in future iterations
      i_to_skip <- c(i_to_skip, row_2_i)
      connected_rows <- cbind(row_1, row_2)  # Combine the original row with the matched PARID row
    }
  }
  
  # Append the connected rows to the new dataset
  STSMs_all_new <- rbind(STSMs_all_new, connected_rows)
}

# Rename the columns of the final dataset to differentiate between the original and matched breakpoints
new_cols <- c(paste0(colnames(stsms_all), 1), paste0(colnames(stsms_all), 2))
colnames(STSMs_all_new) <- new_cols

# Save the final dataset to the specified output path if save_output is set to TRUE
if(save_output == T){
  write_tsv(STSMs_all_new, fs::path(path_results, "Hartwig_all_stsms_info.with_matched_bkpts.no_dup.tsv"))
}


##



