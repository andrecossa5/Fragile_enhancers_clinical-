library(tidyverse)  # Load the tidyverse package for data manipulation
library(VariantAnnotation)  # Load the VariantAnnotation package for handling VCF files
library(GenomicRanges)  # Load the GenomicRanges package for genomic range operations

# Define file paths for input data and results
path_input_main <- fs::path("/hpcnfs/scratch/P_PGP_FRAGILE_ENHANCERS/data")
path_enhancers_ctip <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/downstream/CtIP_enh.hq_signal.clustered.tsv")
path_enhancers_grhl <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/downstream/GRHL_enh.hq_signal.clustered.tsv")
path_results <- fs::path("/hpcnfs/scratch/P_PGP_FRAGILE_ENHANCERS/results/NEW/data/")

# Set a seed for reproducibility
SEED <- 4321
set.seed(SEED)

##

# Define markers and window size
MARKERS <- c("CtIP", "GRHL")
WIN <- 3000

##

## INPUT & PRE-PROCESSING

# Read enhancer files into data frames
suppressMessages({
  enh_ctip <- read_tsv(path_enhancers_ctip) 
  enh_grhl <- read_tsv(path_enhancers_grhl)
  enh_all <- list("CtIP" = enh_ctip, "GRHL" = enh_grhl)
})

# Add chromosome prefix and name column to enhancer data
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

# Convert extended enhancer data frames to GRanges objects
enh_all_ext_gr <- lapply(enh_all_ext, function(df_enh_ext){
  df_enh_ext <- makeGRangesFromDataFrame(df_enh_ext, keep.extra.columns = T)  # Convert to GRanges object
  return(df_enh_ext)
})


##


# Initialize lists to store overlap results for each marker
all.enh.vcf.overlaps <- setNames(vector(mode = "list", length = 2), MARKERS)
all.enh.vcf.overlaps <- list("CtIP" = as.data.frame(matrix(nrow = 0, ncol = 19)), 
                             "GRHL" = as.data.frame(matrix(nrow = 0, ncol = 19)))
cols <- c("seqnames", "start", "end", "width", "strand", "cluster", "summit", "name",
          "seqnames", "start", "end", "width", "strand", "REF", "ALT", "QUAL", "FILTER",  "SAMPLE", "ID")
all.enh.vcf.overlaps <- lapply(all.enh.vcf.overlaps, function(df){colnames(df) <- cols; return(df)})

# Iterate over each directory in the main input path (each directory represents a sample)
for(dir in list.dirs(path_input_main, recursive = F)){
  # Extract the directory name, which represents the sample name
  dir_name <- str_split(dir, pattern = "/", simplify = T)[length(str_split(dir, pattern = "/", simplify = T))]
  cat("\n")
  print(paste0("Iterating over sample: ", dir_name))
  
  # Construct the path to the VCF file for the current sample
  somatic.filename <- paste0(dir, "/purple/", dir_name, ".purple.somatic.vcf.gz")
  # Read the VCF file using the VariantAnnotation package
  somatic.vcf <- readVcf(somatic.filename, genome = "hg19")
  
  # Filter variants to keep only those that pass quality control ("PASS")
  qc_filt <- fixed(somatic.vcf)[, "FILTER"] == "PASS"
  somatic.vcf.filt <- somatic.vcf[qc_filt]
  print("-- Dropping variants NOT passing QC filters --")
  print(paste0("Initial number of variants: ", length(somatic.vcf), " | Number of variants after filtering: ", length(somatic.vcf.filt)))
  
  # Filter to keep only single nucleotide variants (SNVs)
  snv_filt <- (nchar(fixed(somatic.vcf.filt)[, "REF"]) == 1) & unlist((nchar(fixed(somatic.vcf.filt.filt)[, "ALT"]) == 1))
  somatic.vcf.filt.snvs <- somatic.vcf.filt[snv_filt]
  print("-- Retaining only SNVs --")
  print(paste0("Initial number of variants: ", length(somatic.vcf.filt), " | Number of variants after filtering: ", length(somatic.vcf.filt.snvs)))
  
  # Extract the ranges of the filtered SNVs and adjust chromosome names
  somatic.vcf.filt.snvs.gr <- somatic.vcf.filt.snvs@rowRanges
  seqlevels(somatic.vcf.filt.snvs.gr) <- paste0("chr", seqlevels(somatic.vcf.filt.snvs.gr))
  
  # Add metadata to the GRanges object
  df_meta <- cbind(fixed(somatic.vcf.filt.snvs), 
                   data.frame("SAMPLE" = rep(dir_name, length(somatic.vcf.filt.snvs.gr))), 
                   data.frame("PURPLE_AF" = info(somatic.vcf.filt.snvs)[, 'PURPLE_AF']), 
                   data.frame("AF" = geno(somatic.vcf.filt.snvs)$AF[, 2])
  )
  mcols(somatic.vcf.filt.snvs.gr) <- df_meta

  # Compute overlaps between SNVs and enhancers for each marker
  for(marker in MARKERS){
    print(paste0("-- Computing overlaps among SNVs and ", marker, " enhancers --"))
    enh_ext_gr <- enh_all_ext_gr[[marker]]  # Get the GRanges object for the current marker
    
    # Find overlaps between SNVs and enhancer regions
    hits <- findOverlaps(query=enh_ext_gr, subject=somatic.vcf.filt.snvs.gr)
    q <- as.data.frame(enh_ext_gr[queryHits(hits)])  # Extract enhancer data
    s <- cbind(as.data.frame(somatic.vcf.filt.snvs.gr[subjectHits(hits)], row.names=NULL), data.frame("ID" = names(somatic.vcf.filt.snvs.gr[subjectHits(hits)])))  # Extract SNV data
    enh.vcf.overlaps <- cbind(q, s)  # Combine enhancer and SNV data
    
    # Append the overlaps to the list
    all.enh.vcf.overlaps[[marker]] <- rbind(all.enh.vcf.overlaps[[marker]], enh.vcf.overlaps)
  }
}

# Save the overlaps to TSV files for each marker
all.enh.vcf.overlaps[["CtIP"]] %>% write_tsv(., fs::path(path_results, paste0("CtIP_enh.hartwig_snvs.overlap.", WIN, "_3000.tsv")))
all.enh.vcf.overlaps[["GRHL"]] %>% write_tsv(., fs::path(path_results, paste0("GRHL_enh.hartwig_snvs.overlap.", WIN, "_3000.tsv")))
