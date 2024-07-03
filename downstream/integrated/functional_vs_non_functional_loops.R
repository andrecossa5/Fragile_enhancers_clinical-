
#' Script to identify functional and non-functional loops
#' 
#' For now, I am using peaks from a random Chip:  
#' /hpcnfs/data/PGP/madinolfi/Peaks/MCF10A_summit/MCF10A_S23154_NT_190808_hg18_MRE11_peaks.narrowPeak 
#' 
#' @TODO: 
#' - Use peaks from H3k27ac Chip in SCR vs. GRHL2-KD conditions
#' - What if we want to consider as differential only SCR-specific loops?

library(tidyverse)
library(GenomicRanges)
library(rtracklayer)

SEED <- 4321
set.seed(SEED)
pseudocount <- 0.1

path_chip <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/Chip/MCF10A_S23154_NT_190808_hg18_MRE11_peaks.narrowPeak")
path_input_loops <- fs::path("/Users/ieo6983/Desktop/enhancers_project/Analyses/loops/results/2kb/data/tables/2kb_Unified_table.SCR_plus_KD_counts.all_anno_loops.ENH_DEGs_any.tsv")

  
##
  
# Read Chip peaks
peaks_SCR <- rtracklayer::import(path_chip)
# Create fake KD peaks with different signals by changin signalValue column
peaks_KD <- peaks_SCR
peaks_KD$signalValue[1:floor(length(peaks_KD$signalValue)/3)] <- peaks_KD$signalValue[1:floor(length(peaks_KD$signalValue)/3)]+10
m <- min(peaks_KD$signalValue[floor(length(peaks_KD$signalValue)/3):(floor(length(peaks_KD$signalValue)/3)*2)])
peaks_KD$signalValue[floor(length(peaks_KD$signalValue)/3):(floor(length(peaks_KD$signalValue)/3)*2)] <- peaks_KD$signalValue[floor(length(peaks_KD$signalValue)/3):(floor(length(peaks_KD$signalValue)/3)*2)]-m

#' @TODO: 
#' In reality I will have 2 different sets of peaks, called in SCR and KD. So I should:
#' - Obtain the union of the two sets of peaks (some of them will be present in both SCR and KD, other only in 1 condition)
#' - For those peaks present only in one condition, add the peak coordinate in the other condition as well, setting the signalValue to 0
#' - For those peaks overlapping across conditions, extend peaks coordinates to match coordinates of union peaks and keep signalValue of each condition
#' - Compute the FC and log2FC of signalValue in KD vs. SCR for all peaks


# Read loops
# Table contains ALL loops SCR & KD, overlapping a GRHL2 enhancer or a DEG promoter
loops <- read_tsv(path_input_loops)
# Set NA counts to 0, to be able to compute FC
loops$counts.kd[is.na(loops$counts.kd)] <- 0
loops$counts.scr[is.na(loops$counts.scr)] <- 0


##

# Defining differential peaks across conditions - loop FC 
peaks_fc <- peaks_KD$signalValue / peaks_SCR$signalValue
peaks_log2fc <- log2(peaks_fc+pseudocount)
hist(peaks_log2fc, main = "log2(FC) of peaks signal KD vs. SCR")

peaks_union <- data.frame(peaks_SCR[,0])[,1:3]
peaks_union$peaks_log2FC <- peaks_log2fc


##


# Defining differential loops across conditions - loop FC 
loops_fc <- loops$counts.kd / loops$counts.scr
loops_log2fc <- log2(loops_fc+pseudocount)

thresh <- 1
hist(loops_log2fc, main = "log2(FC) of loops counts KD vs. SCR")
abline(v = +thresh, col = "red", lwd = 2, lty = 2)
abline(v = -thresh, col = "red", lwd = 2, lty = 2)

loops$loops_log2FC <- loops_log2fc
#loops$diff <- F
#loops$diff[log2fc > thresh | log2fc < -thresh] <- T


##


# Find overlaps among loop bins and peaks
bin1_gr <- makeGRangesFromDataFrame(loops[, 1:3], seqnames.field = "seqnames1", start.field = "start1", end.field = "end1")
bin2_gr <- makeGRangesFromDataFrame(loops[, 4:6], seqnames.field = "seqnames2", start.field = "start2", end.field = "end2")

peaks_union_gr <- makeGRangesFromDataFrame(peaks_union, keep.extra.columns = T)

findOverlaps(bin1_gr, peaks_union_gr)


