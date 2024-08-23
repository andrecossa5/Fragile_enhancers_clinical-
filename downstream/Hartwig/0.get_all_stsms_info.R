library(tidyverse)  # Load the tidyverse package for data manipulation
library(VariantAnnotation)  # Load the VariantAnnotation package for handling VCF files
library(GenomicRanges)  # Load the GenomicRanges package for genomic range operations

# Define file paths for input and output data
path_input_main <- fs::path("/hpcnfs/scratch/P_PGP_FRAGILE_ENHANCERS/data")
path_results <- fs::path("/hpcnfs/scratch/P_PGP_FRAGILE_ENHANCERS/results/data/")

# Set a seed for reproducibility
SEED <- 4321
set.seed(SEED)


##

# Flag to determine if allele frequency (AF) information should be included
AF <- F

##


# Initialize an empty data frame to store all VCF entries
all.vcf.entries <- as.data.frame(matrix(nrow = 0, ncol = 13))
# Define column names for the VCF entries data frame
cols <- c("seqnames", "start", "end", "width", "strand", "REF", "ALT", "QUAL", "FILTER",  "SAMPLE", "TYPE", "PARID", "ID")
colnames(all.vcf.entries) <- cols

# Iterate over each directory in the main input path (which contains a sample's VCF file)
for(dir in list.dirs(path_input_main, recursive = F)){
  # Extract the directory name, which represents the sample name
  dir_name <- str_split(dir, pattern = "/", simplify = T)[length(str_split(dir, pattern = "/", simplify = T))]
  cat("\n")
  print(paste0("Iterating over sample: ", dir_name))
  
  # Construct the path to the VCF file for the current sample
  somatic.filename <- paste0(dir, "/purple/", dir_name, ".purple.sv.vcf.gz")
  # Read the VCF file using the VariantAnnotation package
  somatic.vcf <- readVcf(somatic.filename, genome = "hg19")
  
  # Filter variants to keep only those that pass quality control ("PASS")
  qc_filt <- fixed(somatic.vcf)[, "FILTER"] == "PASS"
  somatic.vcf.filt <- somatic.vcf[qc_filt]
  print("-- Dropping variants NOT passing QC filters --")
  print(paste0("Initial number of variants: ", length(somatic.vcf), " | Number of variants after filtering: ", length(somatic.vcf.filt)))
  
  # Extract the ranges of the filtered variants and adjust chromosome names to include 'chr'
  somatic.vcf.filt.gr <- somatic.vcf.filt@rowRanges
  seqlevels(somatic.vcf.filt.gr) <- paste0("chr", seqlevels(somatic.vcf.filt.gr))
  
  # Determine the column name for breakpoints Partner ID (PARID)
  if("PARID" %in% colnames(info(somatic.vcf.filt))){
    partner_id <- "PARID"
    partner_id_df <- data.frame("PARID" = info(somatic.vcf.filt)[, partner_id])
  } else {
    partner_id <- "MATEID"
    # Handle cases where MATEID might be empty
    partner_id_df <- data.frame("PARID" = sapply(info(somatic.vcf.filt)[, partner_id], function(x) if (length(x) == 0) NA else unlist(x)))
  }
  
  # Create a metadata data frame with additional information about each variant
  df_meta <- cbind(fixed(somatic.vcf.filt), 
                   data.frame("SAMPLE" = rep(dir_name, length(somatic.vcf.filt.gr)), 
                              "TYPE" = info(somatic.vcf.filt)[, 'EVENTTYPE']))
  
  df_meta <- cbind(df_meta, partner_id_df)
  df_meta$ID <- names(somatic.vcf.filt.gr) # Add an ID column with the variant names
  
  # If AF flag is TRUE, include allele frequency (AF) information
  # For some reason, AF has 2 columns
  if(AF == T){
    df_AF <- lapply(info(somatic.vcf.filt)[, 'PURPLE_AF'], function(x) as.data.frame(t(x))) 
    df_AF <- do.call(rbind, df_AF)
    colnames(df_AF) <- c("PURPLE_AF_1", "PURPLE_AF_2")
    df_meta <- cbind(df_meta, df_AF)
  }
  # Add metadata columns to the genomic ranges object
  mcols(somatic.vcf.filt.gr) <- df_meta
  
  # Convert the genomic ranges object to a data frame and append it to the all.vcf.entries data frame
  somatic.vcf.filt.df <- data.frame(somatic.vcf.filt.gr)
  all.vcf.entries <- rbind(all.vcf.entries, somatic.vcf.filt.df)
}

# Save the combined VCF entries to a TSV file in the specified output path
all.vcf.entries %>% write_tsv(., fs::path(path_results, paste0("Hartwig_all_stsms_info.tsv")))
