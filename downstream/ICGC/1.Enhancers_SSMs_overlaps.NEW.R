
library(tidyverse)
library(GenomicRanges)
library(reshape2)

# Set seed for reproducibility
SEED <- 4321
set.seed(SEED)

WIN <- 3000 # Define window size around enhancer summits
MARKERS <- c("CtIP", "GRHL") # Define markers to analyze 
save_table_overlaps <- F # Option to save overlap tables
location <- "local" # 'local' or 'hpc'

# Load custom genomic functions
source("/Users/ieo6983/Desktop/fragile_enhancer_clinical/utils/functions_genomics.R")

# Define file paths for local environment
path_SSMs <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/genomics/pre_processed_ICGC/simple_somatic_mutation.open.matching_calls.with_AFs.tsv")
path_enhancers_ctip <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/downstream/CtIP_enh.hq_signal.clustered.tsv")
path_enhancers_grhl <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/downstream/GRHL_enh.hq_signal.clustered.tsv")

path_results_data <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/results/ICGC/NEW/enhancers_SSMs_overlaps/data/")
path_results_plots <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/results/ICGC/NEW/enhancers_SSMs_overlaps/plots/")

# Create directories for results if they don't exist
if(!dir.exists(path_results_data)){dir.create(path_results_data, recursive = T)}
if(!dir.exists(path_results_plots)){dir.create(path_results_plots, recursive = T)}

# Override paths and source functions if using HPC environment
if(location == "hpc"){
  source("/hpcnfs/scratch/PGP/Ciacci_et_al/fragile_enhancer_clinical/utils/functions_genomics.R")
  
  path_SSMs <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/data/genomics/pre_processed_ICGC/simple_somatic_mutation.open.matching_calls.with_AFs.tsv")
  path_enhancers_ctip <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/downstream/CtIP_enh.hq_signal.clustered.tsv")
  path_enhancers_grhl <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/downstream/GRHL_enh.hq_signal.clustered.tsv")
  
  path_results_data <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/results/ICGC/NEW/enhancers_SSMs_overlaps/data/")
  path_results_plots <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/results/ICGC/NEW/enhancers_SSMs_overlaps/plots/")
  if(!dir.exists(path_results_data)){dir.create(path_results_data, recursive = T)}
  if(!dir.exists(path_results_plots)){dir.create(path_results_plots, recursive = T)}
  
}


##


## Data Loading ##

# Load SSMs (Simple Somatic Mutations) data
SSMs <- read_tsv(path_SSMs)

# Load enhancer data for both CtIP and GRHL markers
enh_ctip <- read_tsv(path_enhancers_ctip) 
enh_grhl <- read_tsv(path_enhancers_grhl)
enh_all <- list("CtIP" = enh_ctip, "GRHL" = enh_grhl)


## Pre-processing ##

# Standardize enhancer chromosome names and create unique names
enh_all <- lapply(enh_all, function(df){
  df$chrom <- paste0("chr", df$chrom)  # Prefix chromosome with 'chr'
  df$name <- str_c(df$chrom, df$summit, sep=":")  # Create unique name using chrom and summit
  return(df)})

# Filter SSMs to retain only SNVs (Single Nucleotide Variants)
SNVs <- SSMs %>% rename(., c("chromosome_start" = "start", "chromosome_end" = "end")) %>% 
  mutate(., chromosome = paste0("chr", .$chromosome)) %>%  # Prefix chromosome with 'chr'
  dplyr::filter(., end - start == 0)  # Keep only SNVs


##


## Overlap Analysis ##

# Initialize list to store overlaps between enhancers and SNVs
input_variants <- SNVs
output_enh_SNVs <- list()
  
for(marker in MARKERS){
  print(marker)
  mut <- "SNVs"
  
  ## Pre-process enhancer data ##
  marker_enh <- enh_all[[marker]]
  
  # Extend enhancer regions by 'WIN' bp on either side of the summit
  marker_enh$start <- marker_enh$summit - WIN
  marker_enh$end <- marker_enh$summit + WIN
  
  # Convert dataframes to GRanges objects for overlap analysis
  enh_gr <- makeGRangesFromDataFrame(marker_enh, keep.extra.columns = T)
  SNVs_gr <- makeGRangesFromDataFrame(input_variants, keep.extra.columns = T)
  
  
  ##
  
  
  ## Find overlaps between enhancers and SNVs ##
  hits_obj_enh <- findOverlaps(query=enh_gr, subject=SNVs_gr)
  SNVs_sbj <- data.frame(SNVs_gr[subjectHits(hits_obj_enh)])
  colnames(SNVs_sbj)[1:5] <- paste(colnames(SNVs_sbj)[1:5], "sbj", sep = "_")
  
  # Combine enhancer and SNV data based on overlaps
  enh_SNVs <- cbind(data.frame(enh_gr[queryHits(hits_obj_enh)]), SNVs_sbj) 
  print( paste0("Total number of overlaps found: ", dim(enh_SNVs)[1], 
                " for window = ", WIN, " bp") )
  
  # Generate unique ID for each overlap based on SNV position and alleles
  snv_id <- paste(enh_SNVs$seqnames_sbj, enh_SNVs$end_sbj-1, sep=":")
  snv_id <- paste0(snv_id, "_", enh_SNVs$reference_genome_allele, "/", enh_SNVs$mutated_to_allele)
  enh_SNVs$ID <- snv_id
  enh_SNVs <- enh_SNVs %>% relocate(., ID, .before = seqnames_sbj)
  
  # Store overlaps in output list  
  output_enh_SNVs[[marker]] <- enh_SNVs
  
  
  ##
  
  
  ## Optionally save overlaps to file ##
  if(save_table_overlaps == T){
    write_tsv(enh_SNVs, 
              fs::path(path_results_data, paste0("Table_enh_SNVs.", marker, ".all_overlaps.", WIN, "bp_WIN", ".tsv")))
  }
}
  

##

