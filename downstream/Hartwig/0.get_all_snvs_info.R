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


# Initialize an empty data frame to store all VCF entries
all.vcf.entries <- as.data.frame(matrix(nrow = 0, ncol = 13))
# Define column names for the VCF entries data frame
cols <- c("seqnames", "start", "end", "width", "strand", "REF", "ALT", "QUAL", "FILTER",  "SAMPLE", "PURPLE_AF", "AF", "ID")
colnames(all.vcf.entries) <- cols

# Iterate over each directory in the main input path (which contains a sample's VCF file)
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
  
  # Further filter to retain only single nucleotide variants (SNVs)
  snv_filt <- (nchar(fixed(somatic.vcf.filt)[, "REF"]) == 1) & unlist((nchar(fixed(somatic.vcf.filt)[, "ALT"]) == 1))
  somatic.vcf.filt.snvs <- somatic.vcf.filt[snv_filt]
  print("-- Retaining only SNVs --")
  print(paste0("Initial number of variants: ", length(somatic.vcf.filt), " | Number of variants after filtering: ", length(somatic.vcf.filt.snvs)))
  
  # Extract the ranges of the SNVs and adjust chromosome names to include 'chr'
  somatic.vcf.filt.snvs.gr <- somatic.vcf.filt.snvs@rowRanges
  seqlevels(somatic.vcf.filt.snvs.gr) <- paste0("chr", seqlevels(somatic.vcf.filt.snvs.gr))
  
  # Create a metadata data frame with additional information about each variant
  df_meta <- cbind(fixed(somatic.vcf.filt.snvs), 
                   data.frame("SAMPLE" = rep(dir_name, length(somatic.vcf.filt.snvs.gr))), 
                   data.frame("PURPLE_AF" = info(somatic.vcf.filt.snvs)[, 'PURPLE_AF']), 
                   data.frame("AF" = geno(somatic.vcf.filt.snvs)$AF[, 2])
  )
  df_meta$ID <- names(somatic.vcf.filt.snvs.gr)  # Add an ID column with the variant names
  mcols(somatic.vcf.filt.snvs.gr) <- df_meta  # Add metadata columns to the genomic ranges object
  
  # Convert the genomic ranges object to a data frame and append it to the all.vcf.entries data frame
  somatic.vcf.filt.snvs.df <- data.frame(somatic.vcf.filt.snvs.gr)
  all.vcf.entries <- rbind(all.vcf.entries, somatic.vcf.filt.snvs.df)
}


# Save the combined VCF entries to a TSV file in the specified output path
all.vcf.entries %>% write_tsv(., fs::path(path_results, paste0("Hartwig_all_snvs_info.tsv")))
