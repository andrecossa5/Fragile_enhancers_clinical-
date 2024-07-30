
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
library(TFBSTools)
library(BSgenome)
library(JASPAR2020)
library(GenomicRanges)
library(VariantAnnotation)
library(BSgenome.Hsapiens.UCSC.hg19)

SEED <- 4321
set.seed(SEED)
  

dir <- "local" # "local" or "hpc"
if(dir == "hpc"){
  path_enh_SSMs <- list("CtIP" = fs::path(paste0("/hpcnfs/scratch/PGP/Ciacci_et_al/results/ICGC/enhancers_SSMs_overlaps/data/Table_enh_SSMs_CtIP.all_overlaps.3000bp_WIN.tsv")), 
                        "GRHL" = fs::path(paste0("/hpcnfs/scratch/PGP/Ciacci_et_al/results/ICGC/enhancers_SSMs_overlaps/data/Table_enh_SSMs_GRHL.all_overlaps.3000bp_WIN.tsv")))
  path_results <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/results/ICGC/damaging_variants_annotation/")  
} else {
  path_enh_SSMs <- list("CtIP" = fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/results/ICGC/NEW/enhancers_SSMs_overlaps/data/Table_enh_SNVs.CtIP.all_overlaps.3000bp_WIN.tsv"), 
                        "GRHL" = fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/results/ICGC/NEW/enhancers_SSMs_overlaps/data/Table_enh_SNVs.GRHL.all_overlaps.3000bp_WIN.tsv"))
  path_results <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/results/ICGC/NEW/damaging_variants_annotation/")
}

WIN <- 1000
MARKERS <- c("CtIP", "GRHL")
motif_thresh <- 70 # score is typically expressed as a percentage of the maximum possible score for a match to the PWM. 
save_anno <- F
FIMO <- F


##


# GRHL2 motif ID - BaseID: MA1105; MatrixID: MA1105.2 (JASPAR2020 CORE, latest version)
motif_id <- "MA1105.2"
pfm <- getMatrixSet(JASPAR2020, opts = list(ID = motif_id))

# Load variants - enhancers overlap 
enh_SSMs <- list()
enh_gr <- list()
for(marker in MARKERS){
  enh_SSMs[[marker]] <- read_tsv(path_enh_SSMs[[marker]]) %>% suppressMessages()
  
  ## Filter enh_SSMs obj to keep overlaps compatible with WIN
  
  enh_SSMs_filt <- enh_SSMs[[marker]]
  print(paste0("--- Filtering overlaps respecting WIN = ", WIN, " ---"))
  print(paste0("Initial number of overlaps: ", dim(enh_SSMs_filt)[1]))
  
  enh_SSMs_filt$dist <- abs(enh_SSMs_filt$start_sbj - enh_SSMs_filt$summit)
  enh_SSMs_filt <- enh_SSMs_filt %>% dplyr::filter(., dist <= WIN)
  print(paste0("Final number of overlaps after filtering: ", dim(enh_SSMs_filt)[1])); cat("\n")
  
  # Change sequences size as well 
  enh_SSMs_filt$start <- enh_SSMs_filt$start + (3000 - WIN)
  enh_SSMs_filt$end <- enh_SSMs_filt$end - (3000 - WIN)
  enh_SSMs_filt$width <- unique(enh_SSMs_filt$end - enh_SSMs_filt$start)
  enh_SSMs[[marker]] <- enh_SSMs_filt
    
  ##
  
  to_conv <- enh_SSMs[[marker]] %>% dplyr::select(., c(seqnames, start, end, width, summit, name))
  # Convert 50bp extended enhancers into GRanges 
  enh_gr1 <- makeGRangesFromDataFrame(to_conv, keep.extra.columns = T, 
                                      seqnames.field = "seqnames", start.field = "start", end.field = "end", ignore.strand = T)
  
  enh_gr[[marker]] <- enh_gr1
}


##


## Motif scanning across enhancer regions 

# GRHL2 motif has a length of 12 bp, but matchPWM looks for complete matches (of the whole sequence)
pwm <- toPWM(pfm[[1]])

matches_list_markers <- list()
motif_pos_info <- data.frame(matrix(nrow = 1, ncol=5)) 
colnames(motif_pos_info) <- c("enh_name", "motif_id", "motif_start", "motif_end", "motif_score")
motif_pos_info_markers <- list("CtIP" = motif_pos_info, "GRHL" = motif_pos_info)

for(marker in MARKERS){
  
  # Look for GRHL2 motifs matches within enhancer sequences
  seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg19, enh_gr[[marker]])
  min_score_thresh <- paste0(as.character(motif_thresh), "%")
  matches_list <- lapply(seqs, function(seq) {
    matchPWM(pwm@profileMatrix, seq, min.score = min_score_thresh, with.score = TRUE)
  })
  matches_list_markers[[marker]] <- matches_list
  
  # Get motif start-end coordinates in enhancer regions 
  for(i in 1:length(matches_list_markers[[marker]])){
    match <- matches_list_markers[[marker]][[i]]
    if(length(match) == 0){
      # TODO: add enhancer info 
      enh_name <- enh_gr[[marker]][i]$name
      info <- data.frame(matrix(c(enh_name, rep(NA,4)), nrow = 1))
      colnames(info) <- c("enh_name", "motif_id", "motif_start", "motif_end", "motif_score")
      motif_pos_info_markers[[marker]] <- rbind(motif_pos_info_markers[[marker]], info)
    } else {
      # CAREFUL: there can be more than 1 match for each enhancer 
      motif_start <- start(enh_gr[[marker]][i]) + start(match) - 1
      motif_end <-  start(enh_gr[[marker]][i]) + end(match) - 1
      motif_score <- round(mcols(match)$score,2) 
      enh_name <- enh_gr[[marker]][i]$name
      l <- length(motif_start)
      info <- data.frame(matrix(c(rep(enh_name,l), rep(motif_id,l), motif_start, motif_end, motif_score), nrow=l))
      colnames(info) <- c("enh_name", "motif_id", "motif_start", "motif_end", "motif_score")
      motif_pos_info_markers[[marker]] <- rbind(motif_pos_info_markers[[marker]], info)
    }
  }
  
  # Remove duplicated rows - we keep only motif match for each enhancer
  motif_pos_info_markers[[marker]] <- motif_pos_info_markers[[marker]][!duplicated(motif_pos_info_markers[[marker]]), ]
  
  # Print some info
  n_enh_motif <- length(unique(
    motif_pos_info_markers[[marker]] %>% 
      dplyr::filter(., ( !duplicated(.) & !is.na(motif_id) ) ) %>%
      .$enh_name
  ))
  n_enh_tot <- length(unique(enh_SSMs[[marker]]$name))
  print(paste0("Motif match threshold used: ", motif_thresh, " %. For window: ", WIN))
  print(paste0("Fraction of ", marker, " enhancers with a GRHL2 motif within: ", n_enh_motif, " / ", n_enh_tot, 
               " - ~", round(n_enh_motif / n_enh_tot * 100), "%"))
}



##


if(FIMO == T){
  
## FIMO could also be used for motif matching 
# Web tool: https://meme-suite.org/meme/tools/fimo 
# Prepare input for the tool

# Input fasta files 
seqs_ctip <- getSeq(BSgenome.Hsapiens.UCSC.hg19, enh_gr[["CtIP"]])
names(seqs_ctip) <- paste0("sequence", seq_along(seqs_ctip))
seqs_names_conv_ctip <- data.frame("seq_name" = names(seqs_ctip), 
                                   "name" = mcols(enh_gr$CtIP)$name)

seqs_grhl <- getSeq(BSgenome.Hsapiens.UCSC.hg19, enh_gr[["GRHL"]])
names(seqs_grhl) <- paste0("sequence", seq_along(seqs_grhl))
seqs_names_conv_grhl <- data.frame("seq_name" = names(seqs_grhl), 
                                   "name" = mcols(enh_gr$GRHL)$name)

#writeXStringSet(seqs_ctip, filepath = fs::path(path_results, "CtIP_enhancer_regions.WIN_1000.fasta"))
#writeXStringSet(seqs_grhl, filepath = fs::path(path_results, "GRHL_enhancer_regions.WIN_1000.fasta"))

# Motif in MEME format is also needed
# This can be downloaded from: https://jaspar2020.genereg.net/matrix/MA1105.2/ 

# Read FIMO results - table contains duplicated sequences 
seqs_names_conv <- list("CtIP" = seqs_names_conv_ctip, "GRHL" = seqs_names_conv_grhl)
for(marker in MARKERS[2]){
  fimo_results <- read_tsv(fs::path(path_results, paste0("FIMO_Results_info_files.", marker, "_enhancers/fimo.tsv")), comment = "#") %>% suppressMessages()
  fimo_results <- fimo_results[!duplicated(fimo_results),] %>% 
    left_join(., seqs_names_conv[[marker]], by = c("sequence_name" = "seq_name")) %>%
    dplyr::select(., -sequence_name) %>%
    dplyr::filter(., !duplicated(.))
  
  sum(!fimo_results$`p-value` < 0.01) # Only significant matches are returned - in terms of p-value
  sum(fimo_results$`q-value` < 0.05) # However, no significant q-value
  
  # Print some info
  n_enh_tot <- length(unique(enh_SSMs[[marker]]$name))
  n_enh_motif <- length(unique(fimo_results$name))
  cat("\n")
  print("--- Motif scanning results from FIMO ---")
  print(paste0("Fraction of ", marker, " enhancers with a GRHL2 motif within: ", n_enh_motif, " / ", n_enh_tot, " - ~", 
               round(n_enh_motif / n_enh_tot * 100), "%"))
}
}


##


## Annotation of SNVs to damaging levels according to variant position 

enh_SSMs_anno <- list()

for(marker in MARKERS){
  # Prepare SSMs to annotate
  df_enh_var <- enh_SSMs[[marker]]
  df_vars_only <- df_enh_var %>% dplyr::select(., c(ID, seqnames_sbj, start_sbj, end_sbj, icgc_mutation_id))
  df_vars_only_gr <- makeGRangesFromDataFrame(df_vars_only, keep.extra.columns = T, 
                                              seqnames.field = "seqnames_sbj", start.field = "start_sbj", end.field = "end_sbj")

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
  # info already stored in enh_SSMs

  # Annotate
  
  df_enh_var$damaging <- 3
  df_enh_var[hits1, ]$damaging <- 1
  df_enh_var[ (hits2 & !hits1), ]$damaging <- 2
  enh_SSMs_anno[[marker]] <- df_enh_var
  
  # Print some info
  cat("\n")
  print(paste0("Total number of variants for ", marker, " enhancers - WIN = ", WIN, ": ", 
               dim(enh_SSMs_anno[[marker]])[1]))
  t <- table(enh_SSMs_anno[[marker]]$damaging) / sum(table(enh_SSMs_anno[[marker]]$damaging)) * 100
  print("Proportion of damaging vs. non-damaging: ")
  print(paste0( "Damaging level - ", names(t), ": ", round(t, 2), "%" ))
  
  # Save annotated 
  if(save_anno == T){
    enh_SSMs_anno[[marker]] %>% 
      write_tsv(., fs::path(path_results, paste0("Table_enh_SSMs_", marker, ".all_overlaps.", WIN, "bp_WIN.with_damaging_anno.motif_thresh_", motif_thresh, ".tsv")))
  }
}


##

