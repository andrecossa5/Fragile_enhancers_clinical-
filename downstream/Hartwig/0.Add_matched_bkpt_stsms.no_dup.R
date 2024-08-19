
library(tidyverse)

SEED <- 4321
set.seed(SEED)
save_output <- T

path_stsms_all <- fs::path("/hpcnfs/scratch/P_PGP_FRAGILE_ENHANCERS/results/data/Hartwig_all_stsms_info.tsv")
path_results <- fs::path("/hpcnfs/scratch/P_PGP_FRAGILE_ENHANCERS/results/NEW/data/")
  

##


stsms_all <- read_tsv(path_stsms_all)

# Each row contains only 1 bkpt
print("Each row contains only 1 bkpt")
print(sum(!stsms_all$start == stsms_all$end) == 0)

# How many unique SV IDs?
print("Total number of entries - single bkpts: ")
print(dim(stsms_all)[1]) # 871'123

print("Total number of unique IDs and PARIDs:")
print(length(unique(stsms_all$ID))) # 836'226
print(length(unique(stsms_all$PARID))) # 789'061

# Number of IDs is smaller than tot number of entries, and no NAs are present
# Hence, some IDs are duplicated; # While there are NAs among PARIDs
print("Total number of unique IDs and PARIDs:")
sum(is.na(stsms_all$ID))
sum(is.na(stsms_all$PARID))

print("Duplicates in IDs and PARIDs:")
print(sum(duplicated(stsms_all$ID)))
print(sum(duplicated(stsms_all$PARID)))

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


# Add 2nd bkpt to SVS
STSMs_all_new <- data.frame()
i_to_skip <- c()
for(i in 1:dim(stsms_all)[1]){
  print(i)
  
  # Avoid considering rows already checked
  if(i %in% i_to_skip){
    next
  }
  
  # Avoid considering empty rows 
  if( sum(!is.na(stsms_all[i, ])) == 0 ){
    next
  }
  
  row_1 <- stsms_all[i, ] # contains ID and PARID
  
  # PARID might be not indicated
  if(is.na(row_1$PARID)){
    row_2 <- colnames(stsms_all)
    df <- data.frame(matrix(NA, nrow=1, ncol=length(row_2)))
    colnames(df) <- row_2
    connected_rows <- cbind(row_1, df)
    
    STSMs_all_new <- rbind(STSMs_all_new, connected_rows)
  } 
  else {
    # Extract row containing the PARID bkpt 
    # Extract only the PARID row where TYPE and SAMPLE of row_1 and row_2 match
    row_2 <- stsms_all[ (stsms_all$ID == row_1$PARID & stsms_all$TYPE == row_1$TYPE & stsms_all$SAMPLE == row_1$SAMPLE) , ]
    row_2_i <- which(stsms_all$ID == row_1$PARID & stsms_all$TYPE == row_1$TYPE & stsms_all$SAMPLE == row_1$SAMPLE)
    
    # PARID might be not present in the dataset (maybe filtered out)
    if(dim(row_2)[1] == 0){
      row_2 <- colnames(stsms_all)
      df <- data.frame(matrix(NA, nrow=1, ncol=length(row_2)))
      colnames(df) <- row_2
      connected_rows <- cbind(row_1, df)
      
      STSMs_all_new <- rbind(STSMs_all_new, connected_rows)
    }
    # Or, PARID is present in the dataset
    else {
      # Since we are storing row_2 next to row_1 in the new dataset
      # Store index of row_2 and avoid considering it again (creating duplicates)
      i_to_skip <- c(i_to_skip, row_2_i)
      connected_rows <- cbind(row_1, row_2)
    }
  }
  
  STSMs_all_new <- rbind(STSMs_all_new, connected_rows)
}

# Change colnames of final df
new_cols <- c(paste0(colnames(stsms_all), 1), paste0(colnames(stsms_all), 2))
colnames(STSMs_all_new) <- new_cols

if(save_output == T){
  write_tsv(STSMs_all_new, fs::path(path_results, "Hartwig_all_stsms_info.with_matched_bkpts.no_dup.tsv"))
}


##



