library(tidyverse)
library(VariantAnnotation)
library(GenomicRanges)

path_input_main <- fs::path("/hpcnfs/scratch/P_PGP_FRAGILE_ENHANCERS/data")
path_enhancers_ctip <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/downstream/CtIP_enh.hq_signal.clustered.tsv")
path_enhancers_grhl <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/downstream/GRHL_enh.hq_signal.clustered.tsv")
path_results <- fs::path("/hpcnfs/scratch/P_PGP_FRAGILE_ENHANCERS/results/NEW/data/")
  
SEED <- 4321
set.seed(SEED)


##


MARKERS <- c("CtIP", "GRHL")
WIN <- 3000


##

## INPUT & PRE-PROCESSING

# Read input enhancers 
suppressMessages({
  enh_ctip <- read_tsv(path_enhancers_ctip) 
  enh_grhl <- read_tsv(path_enhancers_grhl)
  enh_all <- list("CtIP" = enh_ctip, "GRHL" = enh_grhl)
})

# add summit & name
enh_all <- lapply(enh_all, function(df){
  df$chrom <- paste0("chr", df$chrom)
  df$name <- str_c(df$chrom, df$summit, sep=":")
  return(df)})

# Extend enhancers regions by WIN
enh_all_ext <- lapply(enh_all, function(df_enh_marker){
  df_enh_marker$start <- df_enh_marker$summit - WIN
  df_enh_marker$end <- df_enh_marker$summit + WIN
  return(df_enh_marker)
})

# Convert extended enhancers into GRanges
enh_all_ext_gr <- lapply(enh_all_ext, function(df_enh_ext){
  df_enh_ext <- makeGRangesFromDataFrame(df_enh_ext, keep.extra.columns = T)
  return(df_enh_ext)
})


##


all.enh.vcf.overlaps <- setNames(vector(mode = "list", length = 2), MARKERS)
all.enh.vcf.overlaps <- list("CtIP" = as.data.frame(matrix(nrow = 0, ncol = 19)), 
                             "GRHL" = as.data.frame(matrix(nrow = 0, ncol = 19)))
cols <- c("seqnames", "start", "end", "width", "strand", "cluster", "summit", "name",
          "seqnames", "start", "end", "width", "strand", "REF", "ALT", "QUAL", "FILTER",  "SAMPLE", "ID")
all.enh.vcf.overlaps <- lapply(all.enh.vcf.overlaps, function(df){colnames(df) <- cols; return(df)})

for(dir in list.dirs(path_input_main, recursive = F)){
  dir_name <- str_split(dir, pattern = "/", simplify = T)[length(str_split(dir, pattern = "/", simplify = T))]
  cat("\n")
  print(paste0("Iterating over sample: ", dir_name))
  
  # Read VCF file:
  somatic.filename <- paste0(dir, "/purple/", dir_name, ".purple.somatic.vcf.gz")
  somatic.vcf <- readVcf(somatic.filename, genome = "hg19")
  
  # Keep only variants passing quality filters ("PASS")
  qc_filt <- fixed(somatic.vcf)[, "FILTER"] == "PASS"
  somatic.vcf.filt <- somatic.vcf[qc_filt]
  print("-- Dropping variants NOT passing QC filters --")
  print(paste0("Initial number of variants: ", length(somatic.vcf), " | Number of variants after filtering: ", length(somatic.vcf.filt)))
  
  # Keep only SNVs
  snv_filt <- (nchar(fixed(somatic.vcf.filt)[, "REF"]) == 1) & unlist((nchar(fixed(somatic.vcf.filt)[, "ALT"]) == 1))
  somatic.vcf.filt.snvs <- somatic.vcf.filt[snv_filt]
  print("-- Retaining only SNVs --")
  print(paste0("Initial number of variants: ", length(somatic.vcf.filt), " | Number of variants after filtering: ", length(somatic.vcf.filt.snvs)))
  
  # Extract somatic mutations ranges 
  somatic.vcf.filt.snvs.gr <- somatic.vcf.filt.snvs@rowRanges
  seqlevels(somatic.vcf.filt.snvs.gr) <- paste0("chr", seqlevels(somatic.vcf.filt.snvs.gr))
  
  df_meta <- cbind(fixed(somatic.vcf.filt.snvs), 
                   data.frame("SAMPLE" = rep(dir_name, length(somatic.vcf.filt.snvs.gr))), 
                   data.frame("PURPLE_AF" = info(somatic.vcf.filt.snvs)[, 'PURPLE_AF']), 
                   data.frame("AF" = geno(somatic.vcf.filt.snvs)$AF[, 2])
                   )
  mcols(somatic.vcf.filt.snvs.gr) <- df_meta

  for(marker in MARKERS){
    print(paste0("-- Computing overlaps among SNVs and ", marker, " enhancers --"))
    enh_ext_gr <- enh_all_ext_gr[[marker]]
    
    hits <- findOverlaps(query=enh_ext_gr, subject=somatic.vcf.filt.snvs.gr)
    q <- as.data.frame(enh_ext_gr[queryHits(hits)])
    s <- cbind(as.data.frame(somatic.vcf.filt.snvs.gr[subjectHits(hits)], row.names=NULL),data.frame("ID" = names(somatic.vcf.filt.snvs.gr[subjectHits(hits)])))
    enh.vcf.overlaps <- cbind(q,s)
    
    all.enh.vcf.overlaps[[marker]] <- rbind(all.enh.vcf.overlaps[[marker]], enh.vcf.overlaps)
  }
}

# Save overlaps x marker 
all.enh.vcf.overlaps[["CtIP"]] %>% write_tsv(., fs::path(path_results, paste0("CtIP_enh.hartwig_snvs.overlap.", WIN, "_3000.tsv")))
all.enh.vcf.overlaps[["GRHL"]] %>% write_tsv(., fs::path(path_results, paste0("GRHL_enh.hartwig_snvs.overlap.", WIN, "_3000.tsv")))
