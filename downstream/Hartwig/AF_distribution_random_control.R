

#' Script to add control distribution of AFs from random genomic sequences  
 

cat("\n")
print(" --- Adding a random control to AFs distribution for ICGC and Hartwig --- "); cat("\n")
print("Computing overlap between random genomic regions and ICGC / Hartwig SNVs and saving AFs"); cat("\n")

current_env <- environment()

# Placeholder to save output df with AFs of random regions
ran_AFs_df <- NULL


##


# Create an Isolated Environment
local({
  
suppressMessages({
  
library(tidyverse)
library(GenomicRanges)
source("/hpcnfs/scratch/PGP/Ciacci_et_al/fragile_enhancer_clinical/utils/functions_genomics.R")

SEED <- 4321
set.seed(SEED)

MARKERS <- c("CtIP", "GRHL")
WIN_l <- WIN

path_enhancers_ctip <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/downstream/CtIP_enh.hq_signal.clustered.tsv")
path_enhancers_grhl <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/downstream/GRHL_enh.hq_signal.clustered.tsv")

path_SSMs_icgc <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/data/genomics/pre_processed_ICGC/simple_somatic_mutation.open.matching_calls.with_AFs.tsv")
path_SSMs_hart <- fs::path("/hpcnfs/scratch/P_PGP_FRAGILE_ENHANCERS/results/NEW/data/Hartwig_all_snvs_info.tsv")

path_overlaps_hart <- fs::path("/hpcnfs/scratch/P_PGP_FRAGILE_ENHANCERS/results/NEW/data")
path_overlaps_icgc <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/results/ICGC/enhancers_SSMs_overlaps/NEW/data/")

## Load data

# Read input enhancers 
suppressMessages({
  enh_ctip <- read_tsv(path_enhancers_ctip) 
  enh_grhl <- read_tsv(path_enhancers_grhl)
  enh_all <- list("CtIP" = enh_ctip, "GRHL" = enh_grhl)
})

# Read SSMs
SSMs_icgc <- read_tsv(path_SSMs_icgc)
SSMs_icgc <- SSMs_icgc %>% rename(., c("chromosome_start" = "start", "chromosome_end" = "end")) %>% 
  mutate(., chromosome = paste0("chr", .$chromosome)) %>%
  dplyr::filter(., end - start == 0)
SSMs_hart <- read_tsv(path_SSMs_hart)


AF_ran_all_dfs <- list()
for(marker in MARKERS){
  # Tot. number of enhancers x factor that were overlaped with ICGC and Hartwig somatic SNVs
  n_enh <- dim(enh_all[[marker]])[1]
  
  # Generate same number of random sequences 
  ran_seqs <- generate_random_seqs(SEED, n_enh, win = WIN_l)
  ran_seqs_gr = makeGRangesFromDataFrame(ran_seqs, keep.extra.columns = T)
  
  SSMs_icgc_gr <- makeGRangesFromDataFrame(SSMs_icgc, keep.extra.columns = T)
  SSMs_hart_gr <- makeGRangesFromDataFrame(SSMs_hart, keep.extra.columns = T)
  
  
  ## Compute overlaps 
  
  # Overlap with ICGC
  hits_obj_ran <- findOverlaps(query=ran_seqs_gr, subject=SSMs_icgc_gr)
  SSMs_sbj <- data.frame(SSMs_icgc_gr[subjectHits(hits_obj_ran)])
  colnames(SSMs_sbj)[1:5] <- paste(colnames(SSMs_sbj)[1:5], "sbj", sep = "_")
  ran_seqs_SSMs_icgc <- cbind(data.frame(ran_seqs_gr[queryHits(hits_obj_ran)], SSMs_sbj))

  # Overlap with Hart 
  hits_obj_ran <- findOverlaps(query=ran_seqs_gr, subject=SSMs_hart_gr)
  SSMs_sbj <- data.frame(SSMs_hart_gr[subjectHits(hits_obj_ran)])
  colnames(SSMs_sbj)[1:5] <- paste(colnames(SSMs_sbj)[1:5], "sbj", sep = "_")
  ran_seqs_SSMs_hart <- cbind(data.frame(ran_seqs_gr[queryHits(hits_obj_ran)], SSMs_sbj))
  
  
  ## Prepare dfs
  
  # ICGC
  df_icgc <- ran_seqs_SSMs_icgc %>% dplyr::select(., name, icgc_mutation_id, AF)
  #df_icgc$dist <- NA; df_icgc$cluster <- NA
  colnames(df_icgc) <- c("name", "ID", "AF")
  df_icgc$source <- "ran_icgc"
  
  # HART
  # TODO: chnage back to normal AF
  #df_hart <- ran_seqs_SSMs_hart %>% dplyr::select(., name, ID, AF)
  df_hart <- ran_seqs_SSMs_hart %>% dplyr::select(., name, ID, PURPLE_AF)
  colnames(df_hart)[3] <- "AF"
  
  df_hart$source <- "ran_hart"
  
  # Merge overlaps with random regions for both datasets
  AF_ran_all <- na.omit(rbind(df_icgc, df_hart)) # Many ICGC samples have no AF info 
  AF_ran_all$source <- factor(AF_ran_all$source, levels = c("ran_icgc", "ran_hart"))
  
  AF_ran_all_dfs[[marker]] <- AF_ran_all
}

# Save generated df 
if (exists("AF_ran_all_dfs")) {
  assign("ran_AFs_df", AF_ran_all_dfs, envir = current_env)
}

})
})


##

