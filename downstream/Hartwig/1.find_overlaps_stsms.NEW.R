library(tidyverse)  # Load the tidyverse package for data manipulation
library(VariantAnnotation)  # Load the VariantAnnotation package for handling VCF files
library(GenomicRanges)  # Load the GenomicRanges package for working with genomic ranges

# Define file paths for input data and results
path_input_main <- fs::path("/hpcnfs/scratch/P_PGP_FRAGILE_ENHANCERS/data")
path_enhancers_ctip <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/downstream/CtIP_enh.hq_signal.clustered.tsv")
path_enhancers_grhl <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/downstream/GRHL_enh.hq_signal.clustered.tsv")
path_results <- fs::path("/hpcnfs/scratch/P_PGP_FRAGILE_ENHANCERS/results/NEW/data/")

# Set a seed for reproducibility
SEED <- 4321
set.seed(SEED)

## Define markers and parameters
MARKERS <- c("CtIP", "GRHL")  # Define the markers of interest
WIN <- 3000  # Define the window size for extending enhancer regions
AF <- FALSE  # Flag indicating whether to include allele frequency (AF) information

## INPUT & PRE-PROCESSING

# Read enhancer data from TSV files
suppressMessages({
  enh_ctip <- read_tsv(path_enhancers_ctip) 
  enh_grhl <- read_tsv(path_enhancers_grhl)
  enh_all <- list("CtIP" = enh_ctip, "GRHL" = enh_grhl)  # Store data in a named list
})

# Add chromosome prefix and create a unique name for each enhancer
enh_all <- lapply(enh_all, function(df){
  df$chrom <- paste0("chr", df$chrom)  # Prefix chromosome names with 'chr'
  df$name <- str_c(df$chrom, df$summit, sep=":")  # Create a 'name' column combining chromosome and summit position
  return(df)
})

# Extend enhancer regions by the specified window size
enh_all_ext <- lapply(enh_all, function(df_enh_marker){
  df_enh_marker$start <- df_enh_marker$summit - WIN  # Calculate start position
  df_enh_marker$end <- df_enh_marker$summit + WIN  # Calculate end position
  return(df_enh_marker)
})

# Convert extended enhancer data frames to GRanges objects for overlap analysis
enh_all_ext_gr <- lapply(enh_all_ext, function(df_enh_ext){
  df_enh_ext <- makeGRangesFromDataFrame(df_enh_ext, keep.extra.columns = TRUE)  # Convert to GRanges object
  return(df_enh_ext)
})



##



## Initialize storage for overlap results
# Initialize data frames to store overlaps for each marker
# TODO: check and change ncol if necessary
all.enh.vcf.overlaps <- list("CtIP" = as.data.frame(matrix(nrow = 0, ncol = 21)), 
                             "GRHL" = as.data.frame(matrix(nrow = 0, ncol = 21)))
# Define column names for the overlap results
cols <- c("seqnames", "start", "end", "width", "strand", "cluster", "summit", "name",
          "seqnames", "start", "end", "width", "strand", "REF", "ALT", "QUAL", "FILTER", "SAMPLE", "TYPE", "PARID", "ID")
all.enh.vcf.overlaps <- lapply(all.enh.vcf.overlaps, function(df){colnames(df) <- cols; return(df)})

# Process each directory in the main input path (each directory represents a sample)
for(dir in list.dirs(path_input_main, recursive = F)){
  # Extract the directory name, which represents the sample name
  dir_name <- str_split(dir, pattern = "/", simplify = TRUE)[length(str_split(dir, pattern = "/", simplify = TRUE))]
  cat("\n")
  print(paste0("Iterating over sample: ", dir_name))
  
  # Construct the path to the VCF file for the current sample
  somatic.filename <- paste0(dir, "/purple/", dir_name, ".purple.sv.vcf.gz")
  # Read the VCF file using the VariantAnnotation package
  somatic.vcf <- readVcf(somatic.filename, genome = "hg19")
  
  # Filter variants to keep only those passing quality control ("PASS")
  qc_filt <- fixed(somatic.vcf)[, "FILTER"] == "PASS"
  somatic.vcf.filt <- somatic.vcf[qc_filt]
  print("-- Dropping variants NOT passing QC filters --")
  print(paste0("Initial number of variants: ", length(somatic.vcf), " | Number of variants after filtering: ", length(somatic.vcf.filt)))
  
  # Extract the ranges of the filtered variants and adjust chromosome names
  somatic.vcf.filt.gr <- somatic.vcf.filt@rowRanges
  seqlevels(somatic.vcf.filt.gr) <- paste0("chr", seqlevels(somatic.vcf.filt.gr))
  
  # Get the column indicating breakpoints Partner ID (PARID or MATEID)
  if("PARID" %in% colnames(info(somatic.vcf.filt))){
    partner_id <- "PARID"
    partner_id_df <- data.frame("PARID" = info(somatic.vcf.filt)[, partner_id])
  } else {
    partner_id <- "MATEID"
    partner_id_df <- data.frame("PARID" = sapply(info(somatic.vcf.filt)[, partner_id], function(x) if (length(x) == 0) NA else unlist(x)))
  }
  
  # Create a metadata data frame with additional information
  df_meta <- cbind(fixed(somatic.vcf.filt), 
                   data.frame("SAMPLE" = rep(dir_name, length(somatic.vcf.filt.gr)), 
                              "TYPE" = info(somatic.vcf.filt)[, 'EVENTTYPE']))
  df_meta <- cbind(df_meta, partner_id_df)
  
  # Include allele frequency (AF) information if the AF flag is set to TRUE
  if(AF){
    df_AF <- lapply(info(somatic.vcf.filt)[, 'PURPLE_AF'], function(x) as.data.frame(t(x))) 
    df_AF <- do.call(rbind, df_AF)
    colnames(df_AF) <- c("PURPLE_AF_1", "PURPLE_AF_2")
    df_meta <- cbind(df_meta, df_AF)
  }
  
  # Add metadata to the GRanges object
  mcols(somatic.vcf.filt.gr) <- df_meta
  
  # Compute overlaps between variants and enhancer regions for each marker
  for(marker in MARKERS){
    print(paste0("-- Computing overlaps among SNVs and ", marker, " enhancers --"))
    enh_ext_gr <- enh_all_ext_gr[[marker]]  # Get the GRanges object for the current marker
    
    # Find overlaps between variants and enhancer regions
    hits <- findOverlaps(query=enh_ext_gr, subject=somatic.vcf.filt.gr)
    q <- as.data.frame(enh_ext_gr[queryHits(hits)])  # Extract enhancer data
    s <- cbind(as.data.frame(somatic.vcf.filt.gr[subjectHits(hits)], row.names=NULL), data.frame("ID" = names(somatic.vcf.filt.gr[subjectHits(hits)])))  # Extract variant data
    enh.vcf.overlaps <- cbind(q, s)  # Combine enhancer and variant data
    
    # Append the overlaps to the list for the current marker
    all.enh.vcf.overlaps[[marker]] <- rbind(all.enh.vcf.overlaps[[marker]], enh.vcf.overlaps)
  }
}

# Save the overlap results to TSV files for each marker
all.enh.vcf.overlaps[["CtIP"]] %>% write_tsv(., fs::path(path_results, paste0("CtIP_enh.hartwig_stsms.overlap.WIN_", WIN, ".tsv")))
all.enh.vcf.overlaps[["GRHL"]] %>% write_tsv(., fs::path(path_results, paste0("GRHL_enh.hartwig_stsms.overlap.WIN_", WIN, ".tsv")))
