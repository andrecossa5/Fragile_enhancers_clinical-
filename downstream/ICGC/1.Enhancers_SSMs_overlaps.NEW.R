
library(tidyverse)
library(GenomicRanges)
library(reshape2)

SEED <- 4321
set.seed(SEED)

source("/Users/ieo6983/Desktop/fragile_enhancer_clinical/utils/functions_genomics.R")

path_SSMs <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/genomics/pre_processed_ICGC/simple_somatic_mutation.open.matching_calls.with_AFs.tsv")
path_enhancers_ctip <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/downstream/CtIP_enh.hq_signal.clustered.tsv")
path_enhancers_grhl <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/downstream/GRHL_enh.hq_signal.clustered.tsv")

path_results_data <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/results/ICGC/NEW/enhancers_SSMs_overlaps/data/")
path_results_plots <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/results/ICGC/NEW/enhancers_SSMs_overlaps/plots/")

if(!dir.exists(path_results_data)){dir.create(path_results_data, recursive = T)}
if(!dir.exists(path_results_plots)){dir.create(path_results_plots, recursive = T)}

WIN <- 3000
MARKERS <- c("CtIP", "GRHL") 
save_table_overlaps <- F


##


## Load input files 
SSMs <- read_tsv(path_SSMs)

enh_ctip <- read_tsv(path_enhancers_ctip) 
enh_grhl <- read_tsv(path_enhancers_grhl)
enh_all <- list("CtIP" = enh_ctip, "GRHL" = enh_grhl)


## Pre-process 

# change enhancers name  
enh_all <- lapply(enh_all, function(df){
  df$chrom <- paste0("chr", df$chrom)
  df$name <- str_c(df$chrom, df$summit, sep=":")
  return(df)})

# Filter SSMs (only SNVs)
SNVs <- SSMs %>% rename(., c("chromosome_start" = "start", "chromosome_end" = "end")) %>% 
  mutate(., chromosome = paste0("chr", .$chromosome)) %>%
  dplyr::filter(., end - start == 0)


##


## Find enhancers / SNVs overlaps

input_variants <- SNVs
output_enh_SNVs <- list()
  
for(marker in MARKERS){
  print(marker)
  mut <- "SNVs"
  
  ## Pre/process input data 
  marker_enh <- enh_all[[marker]]
  
  # Extend regions of 'win' from summit
  marker_enh$start <- marker_enh$summit - WIN
  marker_enh$end <- marker_enh$summit + WIN

  # Create GRanges objects 
  enh_gr <- makeGRangesFromDataFrame(marker_enh, keep.extra.columns = T)
  SNVs_gr <- makeGRangesFromDataFrame(input_variants, keep.extra.columns = T)
  
  
  ##
  
  
  ## Find mutations / enhancers overlaps 
  hits_obj_enh <- findOverlaps(query=enh_gr, subject=SNVs_gr)
  SNVs_sbj <- data.frame(SNVs_gr[subjectHits(hits_obj_enh)])
  colnames(SNVs_sbj)[1:5] <- paste(colnames(SNVs_sbj)[1:5], "sbj", sep = "_")
  enh_SNVs <- cbind(data.frame(enh_gr[queryHits(hits_obj_enh)]), SNVs_sbj) 
  print( paste0("Total number of overlaps found: ", dim(enh_SNVs)[1], 
                " for window = ", WIN, " bp") )
  
  # add ID
  snv_id <- paste(paste(enh_SNVs$seqnames_sbj, enh_SNVs$start_sbj-1, sep=":"), enh_SNVs$end_sbj,sep="-")
  snv_id <- paste0(snv_id, "-", enh_SNVs$reference_genome_allele, "/", enh_SNVs$mutated_to_allele)
  enh_SNVs$snv_id <- snv_id
  enh_SNVs <- enh_SNVs %>% relocate(., snv_id, .before = seqnames_sbj)
  
  output_enh_SNVs[[marker]] <- enh_SNVs
  
  
  ##
  
  
  # Save enh_SNVs object for subsequent analyses
  if(save_table_overlaps == T){
    write_tsv(enh_SNVs, 
              fs::path(path_results_data, paste0("Table_enh_SNVs.", marker, ".all_overlaps.", WIN, "bp_WIN", ".tsv")))
  }
}
  



