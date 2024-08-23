

#' Script to add control distribution of AFs from random genomic sequences  
 

cat("\n")
print(" --- Adding a random control to AFs distribution for ICGC and Hartwig --- "); cat("\n")
print("Computing overlap between random genomic regions and ICGC / Hartwig SNVs and saving AFs"); cat("\n")

current_env <- environment()# Capture the current environment

# Placeholder to save output dataframe with AFs of random regions
ran_AFs_df <- NULL


##


# Create an Isolated Environment
local({
  
suppressMessages({  # Suppress messages to reduce output clutter

library(tidyverse)  # Load the tidyverse library for data manipulation
library(GenomicRanges)  # Load GenomicRanges for handling genomic intervals
source("/hpcnfs/scratch/PGP/Ciacci_et_al/fragile_enhancer_clinical/utils/functions_genomics.R") # Source custom genomics functions

# Set a seed for reproducibility
SEED <- 4321
set.seed(SEED)

# Define markers and window length
MARKERS <- c("CtIP", "GRHL")
WIN_l <- WIN # Assumes WIN is predefined

# Define file paths for enhancer data and somatic mutations
path_enhancers_ctip <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/downstream/CtIP_enh.hq_signal.clustered.tsv")
path_enhancers_grhl <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/downstream/GRHL_enh.hq_signal.clustered.tsv")

path_SSMs_icgc <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/data/genomics/pre_processed_ICGC/simple_somatic_mutation.open.matching_calls.with_AFs.tsv")
path_SSMs_hart <- fs::path("/hpcnfs/scratch/P_PGP_FRAGILE_ENHANCERS/results/NEW/data/Hartwig_all_snvs_info.tsv")

path_overlaps_hart <- fs::path("/hpcnfs/scratch/P_PGP_FRAGILE_ENHANCERS/results/NEW/data")
path_overlaps_icgc <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/results/ICGC/enhancers_SSMs_overlaps/NEW/data/")


## Load data

# Read input enhancer datasets 
suppressMessages({
  enh_ctip <- read_tsv(path_enhancers_ctip)  # Load CtIP enhancer data
  enh_grhl <- read_tsv(path_enhancers_grhl)  # Load GRHL enhancer data
  enh_all <- list("CtIP" = enh_ctip, "GRHL" = enh_grhl)  # Combine enhancer data into a list
})

# Read ICGC somatic variants (SSMs)
SSMs_icgc <- read_tsv(path_SSMs_icgc)  # Load ICGC somatic mutation data
SSMs_icgc <- SSMs_icgc %>% 
  rename(., c("chromosome_start" = "start", "chromosome_end" = "end")) %>%  # Rename columns for consistency
  mutate(., chromosome = paste0("chr", .$chromosome)) %>%  # Prefix chromosome names with "chr"
  dplyr::filter(., end - start == 0)  # Filter for SNVs (where start == end)
SSMs_hart <- read_tsv(path_SSMs_hart)  # Load Hartwig somatic mutation data


AF_ran_all_dfs <- list() # Initialize a list to store results for each marker
for(marker in MARKERS){
  # Total number of enhancers for the current marker that overlapped with ICGC and Hartwig somatic SNVs
  n_enh <- dim(enh_all[[marker]])[1]
  
  # Generate the same number of random genomic sequences 
  ran_seqs <- generate_random_seqs(SEED, n_enh, win = WIN_l) # Generate random sequences
  ran_seqs_gr = makeGRangesFromDataFrame(ran_seqs, keep.extra.columns = T) # Convert to GRanges object
  
  # Convert SSMs data to GRanges objects
  SSMs_icgc_gr <- makeGRangesFromDataFrame(SSMs_icgc, keep.extra.columns = T)
  SSMs_hart_gr <- makeGRangesFromDataFrame(SSMs_hart, keep.extra.columns = T)
  
  
  ## Compute overlaps 
  
  # Find overlaps between random sequences and ICGC somatic mutations
  hits_obj_ran <- findOverlaps(query=ran_seqs_gr, subject=SSMs_icgc_gr)
  SSMs_sbj <- data.frame(SSMs_icgc_gr[subjectHits(hits_obj_ran)]) # Extract overlapping SSMs
  colnames(SSMs_sbj)[1:5] <- paste(colnames(SSMs_sbj)[1:5], "sbj", sep = "_") # Rename columns
  ran_seqs_SSMs_icgc <- cbind(data.frame(ran_seqs_gr[queryHits(hits_obj_ran)], SSMs_sbj)) # Combine data

  # Overlap with Hart 
  hits_obj_ran <- findOverlaps(query=ran_seqs_gr, subject=SSMs_hart_gr)
  SSMs_sbj <- data.frame(SSMs_hart_gr[subjectHits(hits_obj_ran)]) # Extract overlapping SSMs
  colnames(SSMs_sbj)[1:5] <- paste(colnames(SSMs_sbj)[1:5], "sbj", sep = "_") # Rename columns
  ran_seqs_SSMs_hart <- cbind(data.frame(ran_seqs_gr[queryHits(hits_obj_ran)], SSMs_sbj))  # Combine data
  
  
  ## Prepare dataframes
  
  
  # Prepare dataframe for ICGC overlaps
  df_icgc <- ran_seqs_SSMs_icgc %>% dplyr::select(., name, icgc_mutation_id, AF)  # Select relevant columns
  # df_icgc$dist <- NA; df_icgc$cluster <- NA  # Placeholder columns (commented out)
  colnames(df_icgc) <- c("name", "ID", "AF")  # Rename columns
  df_icgc$source <- "ran_icgc"  # Add source label
  
  # Prepare dataframe for Hartwig overlaps
  df_hart <- ran_seqs_SSMs_hart %>% dplyr::select(., name, ID, AF)  # Original column selection (commented out)
  df_hart$source <- "ran_hart"  # Add source label  
  
  # Merge overlaps with random regions for both datasets
  AF_ran_all <- na.omit(rbind(df_icgc, df_hart))  # Remove rows with missing values and combine datasets
  AF_ran_all$source <- factor(AF_ran_all$source, levels = c("ran_icgc", "ran_hart"))  # Set factor levels for source
  
  AF_ran_all_dfs[[marker]] <- AF_ran_all # Store result in list
}

# Save the generated dataframe to the global environment 
if (exists("AF_ran_all_dfs")) {
  assign("ran_AFs_df", AF_ran_all_dfs, envir = current_env)
}

})
})


##

