
library(tidyverse)

# Set a random seed for reproducibility
SEED <- 4321
set.seed(SEED)

set <- "NEW"  # NEW or OLD, for new or old set of enhancers 
save_output <- F # Specify whether to save the output or not
location <- "local" # 'local' or 'hpc'

# Define the path to the loops data based on the location and set type
path_loops <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/results/integrated/", set)

# If the location is 'hpc', update the path to the HPC location
if(location == "hpc"){
  path_loops <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/results/integrated/NEW/")
}
                             

##


# All loops contain at least 1 enhancer #

# Annotate enhancers for all loops containing at least one enhancer
for(kb in c(2,4)){
  
  # Print a message indicating the kb size of loops being processed
  print(paste0("Annotating enhancers from ", kb, " kb loops"))
  
  # Read the loops data from the specified file
  loops_table <-  read_tsv(fs::path(path_loops, 
    sprintf("%skb/data/tables/%skb_Unified_table.SCR_plus_KD_counts.all_anno_loops.ENH_DEGs_any.tsv", kb, kb))) %>% 
    suppressMessages()
  
  # Filter the data to select SCR-specific loops linked to DOWN DEGs
  loi <- loops_table %>% 
    filter(is.na(counts.kd)) %>% #Select SCR-specific loops 
    filter((!is.na(name1) & !is.na(gene_name2) & DE2 == "Down") | 
           (!is.na(name2) & !is.na(gene_name1) & DE1 == "Down")) #Select ENH-DOWN loops only
  
  # Create a unique identifier for each loop
  loi$unique_loop_id <- str_c(loi$seqnames1, loi$start1, loi$end2, sep = "_")
  
  # Collect enhancer names from the filtered loops
  enh_bin1 <- loi[!is.na(loi$gene_name2), ]$name1 
  enh_bin2 <- loi[!is.na(loi$gene_name1), ]$name2
  enh_oi <- c(enh_bin1[!is.na(enh_bin1)], enh_bin2[!is.na(enh_bin2)])
  enh_oi <- data.frame("name" = unique(enh_oi))
  
  # Print the total number of unique enhancers identified
  print(paste0("Total number of enhancers from SCR-specific-Down loops, ", kb, " kb: ", dim(enh_oi)[1]))
  
  # If saving output is enabled, create directories and save the results  
  if(save_output == T){
    path_results <- fs::path(path_loops, sprintf("%skb/data/anno_enhancers/", kb))
    # Create directory if it does not exist
    if(!dir.exists(path_results)){dir.create(path_results)}
    
    # Save enhancers 
    file_name <- sprintf("%skb_GRHL2_enhancers.from_SCR_specific_loops.linked_to_DOWN_DEGs.tsv", kb)
    enh_oi %>% write_tsv(., fs::path(path_results, file_name))
    
    # Save the filtered loops data
    loi %>% write_tsv(., fs::path(path_results, paste0(kb, "kb_functional_loops.aka_SCR_specific_loops_linked_to_DOWN_DEGs.tsv")))
  }
}


##


