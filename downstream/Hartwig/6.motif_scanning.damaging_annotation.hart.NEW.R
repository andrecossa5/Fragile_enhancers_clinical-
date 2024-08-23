
#' This script can be used to classify mutations as damaging AFTER overlap among enhancers and SSMs/SNVs
#' 
#' @description
#' 3 "damaging" levels for mutations are defined: 
#' - 1 = SNV falls within motif;
#' - 2 = SNV falls outside motif, within a +/- 10 bp window around motif;
#' - 3 = not damaging - SNV falls outside window around motif but within +/- WIN from peak summit.
#' 


##


library(tidyverse)
library(TFBSTools) # For working with position weight matrices (PWM)
library(BSgenome) # For genomic sequences
library(JASPAR2020) # For JASPAR database of PWMs
library(GenomicRanges) # For handling genomic intervals
library(VariantAnnotation) # For working with variant data
library(BSgenome.Hsapiens.UCSC.hg19) # Genomic sequences for human hg19

# Set random seed for reproducibility
SEED <- 4321
set.seed(SEED)

# Define file paths for input data and output results
path_enh_SSMs <- list(
  "CtIP" = fs::path(paste0("/hpcnfs/scratch/P_PGP_FRAGILE_ENHANCERS/results/NEW/data/CtIP_enh.hartwig_snvs.overlap.WIN_3000.tsv")), 
  "GRHL" = fs::path(paste0("/hpcnfs/scratch/P_PGP_FRAGILE_ENHANCERS/results/NEW/data/GRHL_enh.hartwig_snvs.overlap.WIN_3000.tsv"))
)
path_results <- fs::path("/hpcnfs/scratch/P_PGP_FRAGILE_ENHANCERS/results/NEW/data/damaging_variants_annotation/")  
if(!dir.exists(path_results)){dir.create(path_results, recursive = T)}

# Define parameters
WIN <- 1000  # Window size around motifs for classification
MARKERS <- c("CtIP", "GRHL")  # List of markers to process
motif_thresh <- 70  # Threshold score percentage for motif matching
save_anno <- T  # Whether to save the annotated results


##


# Define the motif ID for GRHL2 and retrieve the position weight matrix (PWM) from JASPAR2020
# GRHL2 motif ID - BaseID: MA1105; MatrixID: MA1105.2 (JASPAR2020 CORE, latest version)
motif_id <- "MA1105.2"
pfm <- getMatrixSet(JASPAR2020, opts = list(ID = motif_id))

# Initialize lists to store data for enhancers and genomic ranges
enh_SSMs <- list()
enh_gr <- list()

# Process each marker
for(marker in MARKERS){
  # Define column names for reading input files
  column_names <- c(
    paste0(c("seqnames", "start", "end", "width", "strand", "name", "summit", "cluster"), "_enh"), 
    paste0(c("seqnames", "start", "end", "width", "strand", "REF", "ALT", "QUAL", "FILTER", "SAMPLE", "PURPLE_AF", "AF", "ID"), "_snv")
  )
  
  # Read the data from the input file
  enh_SSMs[[marker]] <- read_tsv(path_enh_SSMs[[marker]])
  colnames(enh_SSMs[[marker]]) <- column_names
  
  ## Filter data to respect the window size around motifs
  
  enh_SSMs_filt <- enh_SSMs[[marker]]
  print(paste0("--- Filtering overlaps respecting WIN = ", WIN, " ---"))
  print(paste0("Initial number of overlaps: ", dim(enh_SSMs_filt)[1]))
  
  # Calculate distance between SNV start and enhancer summit
  enh_SSMs_filt$dist <- abs(enh_SSMs_filt$start_snv - enh_SSMs_filt$summit_enh)
  
  # Keep only overlaps within the defined window size
  enh_SSMs_filt <- enh_SSMs_filt %>% dplyr::filter(., dist <= WIN)
  print(paste0("Final number of overlaps after filtering: ", dim(enh_SSMs_filt)[1])); cat("\n")
  
  # Adjust enhancer coordinates to match the window size
  enh_SSMs_filt$start_enh <- enh_SSMs_filt$start_enh + (3000 - WIN)
  enh_SSMs_filt$end_enh <- enh_SSMs_filt$end_enh - (3000 - WIN)
  enh_SSMs_filt$width_enh <- unique(enh_SSMs_filt$end_enh - enh_SSMs_filt$start_enh)
  enh_SSMs[[marker]] <- enh_SSMs_filt
  
  ##
  
  # Convert enhancer data to GRanges format
  to_conv <- enh_SSMs[[marker]] %>% dplyr::select(., c(seqnames_enh, start_enh, end_enh, width_enh, summit_enh, name_enh))
  enh_gr1 <- makeGRangesFromDataFrame(to_conv, keep.extra.columns = T, 
                                      seqnames.field = "seqnames_enh", start.field = "start_enh", end.field = "end_enh", ignore.strand = T)
  
  enh_gr[[marker]] <- enh_gr1
}


##


## Motif scanning

# Convert the PWM to a format suitable for matching
pwm <- toPWM(pfm[[1]])

# Initialize lists to store motif matches and position information
matches_list_markers <- list()
motif_pos_info <- data.frame(matrix(nrow = 1, ncol=5)) 
colnames(motif_pos_info) <- c("enh_name", "motif_id", "motif_start", "motif_end", "motif_score")
motif_pos_info_markers <- list("CtIP" = motif_pos_info, "GRHL" = motif_pos_info)

# Scan for motifs in each marker
for(marker in MARKERS){
  
  # Extract sequences for the enhancer regions
  seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg19, enh_gr[[marker]])
  
  # Match the PWM to the sequences
  min_score_thresh <- paste0(as.character(motif_thresh), "%")
  matches_list <- lapply(seqs, function(seq) {
    matchPWM(pwm@profileMatrix, seq, min.score = min_score_thresh, with.score = TRUE)
  })
  matches_list_markers[[marker]] <- matches_list
  
  # Record motif match information
  for(i in 1:length(matches_list_markers[[marker]])){
    match <- matches_list_markers[[marker]][[i]]
    if(length(match) == 0){
      # If no matches, record NA values
      enh_name <- enh_gr[[marker]][i]$name_enh
      info <- data.frame(matrix(c(enh_name, rep(NA,4)), nrow = 1))
      colnames(info) <- c("enh_name", "motif_id", "motif_start", "motif_end", "motif_score")
      motif_pos_info_markers[[marker]] <- rbind(motif_pos_info_markers[[marker]], info)
    } else {
      # Record start, end, and score for each motif match
      motif_start <- start(enh_gr[[marker]][i]) + start(match) - 1
      motif_end <-  start(enh_gr[[marker]][i]) + end(match) - 1
      motif_score <- round(mcols(match)$score,2) 
      enh_name <- enh_gr[[marker]][i]$name_enh
      l <- length(motif_start)
      info <- data.frame(matrix(c(rep(enh_name,l), rep(motif_id,l), motif_start, motif_end, motif_score), nrow=l))
      colnames(info) <- c("enh_name", "motif_id", "motif_start", "motif_end", "motif_score")
      motif_pos_info_markers[[marker]] <- rbind(motif_pos_info_markers[[marker]], info)
    }
  }
  
  # Remove duplicated rows - we keep only motif match for each enhancer
  motif_pos_info_markers[[marker]] <- motif_pos_info_markers[[marker]][!duplicated(motif_pos_info_markers[[marker]]), ]
  
  # Print information about motif matches
  n_enh_motif <- length(unique(
    motif_pos_info_markers[[marker]] %>% 
      dplyr::filter(., ( !duplicated(.) & !is.na(motif_id) ) ) %>%
      .$enh_name
  ))
  n_enh_tot <- length(unique(enh_SSMs[[marker]]$name_enh))
  print(paste0("Motif match threshold used: ", motif_thresh, " %. For window: ", WIN))
  print(paste0("Fraction of ", marker, " enhancers with a GRHL2 motif within: ", n_enh_motif, " / ", n_enh_tot, 
               " - ~", round(n_enh_motif / n_enh_tot * 100), "%"))
}


##


## Annotation of SNVs to damaging levels according to variant position 

# Initialize list to store annotated data
enh_SSMs_anno <- list()

# Annotate each marker
for(marker in MARKERS){
  # Prepare data for annotation
  df_enh_var <- enh_SSMs[[marker]] 
  df_vars_only <- df_enh_var[, c("seqnames_snv", "start_snv", "end_snv", "ID_snv")] # icgc: seqnames, start, end, icgc_mut_id
  df_vars_only_gr <- makeGRangesFromDataFrame(df_vars_only, keep.extra.columns = T, 
                                              seqnames.field = "seqnames_snv", start.field = "start_snv", end.field = "end_snv")
  
  # Prepare motifs coordinates for annotation
  motifs_coords <- na.omit(motif_pos_info_markers[[marker]])
  motifs_coords$motif_end <- as.numeric(motifs_coords$motif_end); motifs_coords$motif_start <- as.numeric(motifs_coords$motif_start)
  motifs_coords$chrom <- str_split(motifs_coords$enh_name, ":", simplify = T)[,1]
  motifs_coords_gr <- makeGRangesFromDataFrame(motifs_coords, keep.extra.columns = T)
  
  ## Intersect
  
  # Level 1 - SNV falls WITHIN motif sequence
  hits1 <- !is.na(findOverlaps(df_vars_only_gr, motifs_coords_gr, select = "arbitrary")) # if hits1 == TRUE, SNV falls within motif seq
  
  # Level 2 - SNV falls within motif sequence +/- 10 bp
  motifs_coords_gr_ext <- motifs_coords_gr
  start(motifs_coords_gr_ext) <- start(motifs_coords_gr_ext) - 10
  end(motifs_coords_gr_ext) <- end(motifs_coords_gr_ext) + 10
  hits2 <- !is.na(findOverlaps(df_vars_only_gr, motifs_coords_gr_ext, select = "arbitrary")) # if hits1 == TRUE, SNV falls within motif seq +/- 10 bp 
  
  # Level 3 - SNV does NOT fall within or around motif sequences, but within +/- 50 bp from peak summit
  # Information already stored in enh_SSMs
  
  # Annotate
  
  # Annotate variants based on overlap levels
  df_enh_var$damaging <- 3  # Default level
  df_enh_var[hits1, ]$damaging <- 1  # Level 1: Within motif
  df_enh_var[ (hits2 & !hits1), ]$damaging <- 2  # Level 2: Within extended region
  enh_SSMs_anno[[marker]] <- df_enh_var
  
  # Print information about annotation
  cat("\n")
  print(paste0("Total number of variants for ", marker, " enhancers - WIN = ", WIN, ": ", 
               dim(enh_SSMs_anno[[marker]])[1]))
  t <- table(enh_SSMs_anno[[marker]]$damaging) / sum(table(enh_SSMs_anno[[marker]]$damaging)) * 100
  print("Proportion of damaging vs. non-damaging: ")
  print(paste0( "Damaging level - ", names(t), ": ", round(t, 2), "%" ))
  
  # Save the annotated data if specified
  if(save_anno == T){
  enh_SSMs_anno[[marker]] %>% 
    write_tsv(., fs::path(path_results, paste0("Table_enh_SSMs_", marker, ".all_overlaps.", WIN, "bp_WIN.with_damaging_anno.motif_thresh_", motif_thresh, ".tsv")))
}
}


##



