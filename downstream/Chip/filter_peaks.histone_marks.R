
library(tidyverse)
library(rowr)
library(GenomicRanges)


##


#' @CHIPs inventory
#' 
#' @CtIP
#' _Keep_:
#' - MCF10A_S30627_NT_201120_hg19_CtIP_Total.bigWig
#' - MCF10A_S32649_NT_210226_hg19_CtIP.bigWig
#' - MCF10A_S40470_EMTnegCTRL2_220211_hg19_CtIP_Tot.bigWig
#' - MCF10A_S43554_SCR_220801_hg19_CtIP_Tot.bigWig
#' 
#' _Remove_:
#' - All S327P: --> different signal
#'   MCF10A_S31641_NT_210114_hg19_CtIP_S327P.bigWig, MCF10A_S43555_SCR_220801_hg19_CtIP_S327P.bigWig, MCF10A_S43555_SCR_220801_hg19_CtIP_S327P.bigWig
#' - S327 (without P): --> different signal
#'   MCF10A_S32696_NT_210226_hg19_CtIP_S327.bigWig  
#' - MCF10A_S30627_NT_201112_hg19_CtIP_Total.bigWig --> bad quality 
#' 
#' @GRHL1/2
#' _Remove_:
#' - MCF10A_S40468_SCR_220211_hg19_GRHL1.bigWig --> Bad signal 
#' 
#' _Keep_:
#' - MCF10A_S34023_NT_210512_hg19_GRHL2.bigWig
#' - MCF10A_S36604_NT_210818_hg19_GRHL1.bigWig


##


SEED <- 4321
set.seed(SEED)

MARKERS <- c("CtIP", "GRHL")

# Paths CtIP/GRHL peaks
path_input_peaks <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05//macs2/")
path_peaks_merged <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05//downstream/peaks_union/")

# Paths histone marks 
path_peaks_k4me3 <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/Chip/Chip_for_clusters/results/HISTONE_MARKS/macs2/MCF10A_S48547_230308_hg19_H3K4me3_peaks.narrowPeak")
path_peaks_k4me1 <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/Chip/Chip_for_clusters/results/HISTONE_MARKS/macs2/MCF10A_S47229_SCR_221230_hg19_H3K4me1_peaks.narrowPeak")
path_peaks_k27ac <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/Chip/Chip_for_clusters/results/HISTONE_MARKS/macs2/MCF10A_S43557_SCR_220801_hg19_H3K27ac_peaks.narrowPeak")

# Paths output 
path_output_k4me3 <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/Chip/Chip_for_clusters/results/HISTONE_MARKS/macs2/")
path_output_merged <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05//downstream/peaks_union/")

location <- "local"
if(location == "hpc"){
  path_input_peaks <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/bwa/mergedLibrary/macs2/narrowPeak/")  
  path_output_merged <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/peaks_union/")
}


# metric1 and metric2 could be pValue and qValue
narrowPeak_colnames <- c("chrom", "start", "end", "name", "score", "strand", 
                         "signal_value", "metric1", "metric2", "summit")
narrowPeak_colnames_nostrand <- c("chrom", "start", "end", "name", "score", 
                                  "signal_value", "metric1", "metric2", "summit")


##


## Plot tot. number of peaks x Chip + Peaks annotations x Chip 

# Total number of peaks
num_peaks <- read_tsv(fs::path(path_input_peaks, "qc/macs2_peak.summary.txt")) %>% suppressMessages()
num_peaks[num_peaks$measure == "fold", ] %>% 
  ggplot(., aes(x = sample, y = num_peaks, fill = sample))+
  geom_bar(stat = "identity")+
  theme_light()+
  labs(title = "Number of peaks", y = "n. peaks")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "none")

# Peaks annotations 
anno <- read_tsv(fs::path(path_input_peaks, "qc/macs2_annotatePeaks.summary.txt")) %>% suppressMessages()
anno[is.na(anno)] <- 0
anno <- anno %>% mutate(., tot = Intergenic+TTS+exon+intron+`promoter-TSS`)

anno %>% mutate(., across(Intergenic:`promoter-TSS`, ~ . / tot)) %>% 
  dplyr::select(., -tot) %>%
  pivot_longer(-sample) %>% 
  ggplot(., aes(x=sample, y = value, fill = name))+
  geom_bar(stat = "identity", position = "stack")+
  theme_light()+
  labs(title = "Peaks annotation", y = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))


##


# Filter out potential enhancer peaks according to annotation
anno_of_int <- tolower(c("Intergenic", "intron")) # leaving out: TSS, Promoter, exon

peaks_all <- list.files(path_input_peaks)[endsWith(list.files(path_input_peaks), "peaks.narrowPeak")]
peaks_anno_all <- list.files(path_input_peaks)[endsWith(list.files(path_input_peaks), "annotatePeaks.txt")]

for(marker in MARKERS){
  cat("\n")
  print(marker)
  
  peaks_anno_marker <- grep(marker, peaks_anno_all, value = T)
  print(peaks_anno_marker)
  
  for(anno_file in peaks_anno_marker){
    sample_name <- str_split(anno_file, "\\.")[[1]][1]
    
    peaks_anno <- read_tsv(fs::path(path_input_peaks, anno_file)) %>% suppressMessages()
    colnames(peaks_anno)[1] <- "PeakID"
    peaks_anno$simple_annotation <- str_split(peaks_anno$Annotation, pattern = c("-"), simplify = T)[,1] %>% str_split(., pattern=" ", simplify = T) %>% .[,1]
    
    # Select only peaks that could represent enhancers and filter bed file 
    peaks_anno_filt <- peaks_anno %>% dplyr::filter(., tolower(simple_annotation) %in% anno_of_int)
    
    peaks_bed <- read_tsv(fs::path(path_input_peaks, peaks_all[startsWith(peaks_all, sample_name)]), 
                          col_names = narrowPeak_colnames) %>% suppressMessages() 
    peaks_bed_filt <- peaks_bed[peaks_bed$name %in% peaks_anno_filt$PeakID, ]
    print(dim(peaks_bed_filt)[1] == dim(peaks_anno_filt)[1])
    
    #write_tsv(peaks_bed_filt, fs::path(path_input_peaks, paste0(sample_name, ".filt.narrowPeak" )), col_names = F)
  }
}


##


# --> Merge peaks with bedtools 
# TODO: add path to script


##


## Check distribution of H3K4me3 signal 

peaks_k4me3 <- read_tsv(path_peaks_k4me3, col_names = narrowPeak_colnames)

# signal_value (enrichment)
ggplot(peaks_k4me3, aes(x = signal_value))+
  geom_density(adjust = 0.5)+
  geom_vline(xintercept = 40, color ="red")
summary(peaks_k4me3$signal_value)

# score
ggplot(peaks_k4me3, aes(x = score))+
  geom_density(adjust = 0.5)+
  geom_vline(xintercept = median(peaks_k4me3$score), color ="red")+
  labs(title="H3K4me3 score distribution")
summary(peaks_k4me3$score)

# filter peaks with a score > median(score): these will be considered as 'actual' H3K4me3 peaks 
ms <- median(peaks_k4me3$score)
peaks_k4me3_filt <- peaks_k4me3 %>% dplyr::filter(., score >= ms)

path_output_k4me3_filt <- fs::path(path_output_k4me3, "MCF10A_S48547_230308_hg19_H3K4me3_peaks.filt_score_above_median.narrowPeak")
#peaks_k4me3_filt %>% write_tsv(., path_output_k4me3_filt, col_names = FALSE)


##


## Check distribution of H3K4me1 signal 

peaks_k4me1 <- read_tsv(path_peaks_k4me1, col_names = narrowPeak_colnames)

# signal_value (enrichment)
ggplot(peaks_k4me1, aes(x = signal_value))+
  geom_density(adjust = 1)
summary(peaks_k4me1$signal_value)

# score
ggplot(peaks_k4me1, aes(x = score))+
  geom_density(adjust = 1)+
  geom_vline(xintercept = median(peaks_k4me1$score), color ="red")+
  labs(title = "H3K4me1 score distribution")
summary(peaks_k4me1$score)


##


## Some stats for histone marks

suppressMessages({
  peaks_k27ac <- read_tsv(path_peaks_k27ac, col_names = narrowPeak_colnames)
  ctip_merged <- read_tsv(fs::path(path_peaks_merged, paste0("/merged_peaks.hg19_CtIP.hq_signal.narrowPeak")), col_names = narrowPeak_colnames)
  grhl_merged <- read_tsv(fs::path(path_peaks_merged, paste0("/merged_peaks.hg19_GRHL.hq_signal.narrowPeak")), col_names = narrowPeak_colnames)
})

# Tot. number of peaks
peaks_dims <- data.frame(
  "factor" = c("ctip_merged", "grhl_merged", "H3K27ac", "H3K4me1", "H3K4me3"), 
  "n_peaks" = c(dim(ctip_merged)[1], dim(grhl_merged)[1], dim(peaks_k27ac)[1], dim(peaks_k4me1)[1], dim(peaks_k4me3)[1])
)

peaks_dims %>% 
  ggplot(., aes(x = factor, y = n_peaks, fill = factor))+
  geom_bar(stat = "identity")+
  labs(title = "Number of peaks", x = "")


# Peaks lengths
peaks_lengths <- data.frame(rowr::cbind.fill(ctip_merged$end - ctip_merged$start, 
                 grhl_merged$end - grhl_merged$start, 
                 peaks_k27ac$end - peaks_k27ac$start,
                 peaks_k4me1$end - peaks_k4me1$start, 
                 peaks_k4me3$end - peaks_k4me3$start))
colnames(peaks_lengths) <- peaks_dims$factor

peaks_lengths %>% 
  pivot_longer(everything(), names_to = "factor", values_to = "value") %>%
  ggplot(., aes(x = factor,y = value, fill = factor))+
  geom_boxplot()+
  theme_light()+
  labs(title = "Peaks length distribution", x = "", y = "bp")

peaks_lengths %>% 
  pivot_longer(everything(), names_to = "factor", values_to = "value") %>%
  ggplot(., aes(x = factor,y = value, fill = factor))+
  scale_y_log10()+
  geom_boxplot()+
  theme_light()+
  labs(title = "Peaks length distribution", x = "", y = "log10(bp)")

print("Peaks lengths stats for CtIP:")
summary(peaks_lengths$ctip_merged)
print("Peaks lengths stats for GRHL:")
summary(peaks_lengths$grhl_merged)


##


## Filter CtIP & GRHL1/2 merged based on :
# - Overlap with H3k27ac peaks and H3K4me1 peaks (marker for active enhancers)
# - Absence of overlap with H3K4me3 peaks (marker for active promoters)

suppressMessages({
  peaks_k4me1 <- read_tsv(path_peaks_k4me1, col_names = narrowPeak_colnames) %>%
    makeGRangesFromDataFrame(., keep.extra.columns=T)
  peaks_k4me3_filt <- read_tsv(path_output_k4me3_filt, col_names = narrowPeak_colnames) %>% 
    makeGRangesFromDataFrame(., keep.extra.columns=T)
  peaks_k27ac <- read_tsv(path_peaks_k27ac, col_names = narrowPeak_colnames) %>%
    makeGRangesFromDataFrame(., keep.extra.columns=T)
})

for(marker in MARKERS){
  cat("\n")
  peaks_union <- read_tsv(fs::path(path_peaks_merged, paste0("/merged_peaks.hg19_", marker, ".hq_signal.narrowPeak")), 
                          col_names = narrowPeak_colnames) %>% suppressMessages()
  peaks_union$strand <- NULL # remove since it can me ambiguous  
  peaks_union <- makeGRangesFromDataFrame(peaks_union, keep.extra.columns=T)
  print(paste0("Total n. of ", marker, " peaks: ", length(peaks_union)))
  
  # Remove peaks overlapping with k4me3 peaks
  filt_k4me3 <- is.na(findOverlaps(peaks_union, peaks_k4me3_filt, select = "arbitrary"))
  peaks_union_filt1 <- peaks_union[filt_k4me3] 
  print(paste0("Total n. of ", marker, " peaks - after k4me3 filtering: ", length(peaks_union_filt1)))
  
  # Keep only peaks overlapping with H3k27ac & H3k4me1 peaks
  filt_k27ac <- !is.na(findOverlaps(peaks_union_filt1, peaks_k27ac, select = "arbitrary"))
  filt_k4me1 <- !is.na(findOverlaps(peaks_union_filt1, peaks_k4me1, select = "arbitrary")) 

  peaks_union_filt2 <- peaks_union_filt1[filt_k27ac & filt_k4me1]
  print(paste0("Total n. of ", marker, " peaks - after k4me1 & k27ac filtering: ", length(peaks_union_filt2)))
  
  # Save new set of filtered peaks 
  peaks_union_filt2_df <- data.frame(peaks_union_filt2) %>% dplyr::select(., -c(4,5))
  #peaks_union_filt2_df %>% write_tsv(., fs::path(path_output_merged, paste0("merged_peaks.hg19_", marker, ".filtered_K27ac_k4me1_k4me3", ".hq_signal.narrowPeak")),col_names = F)
}


##


## Add new summit 

# Option 1: add summit as the middle point of each peak 
for(marker in MARKERS){
  peaks_union <- read_tsv(fs::path(path_peaks_merged, paste0("merged_peaks.hg19_", marker, ".filtered_K27ac_k4me1_k4me3.hq_signal.narrowPeak")), 
                          col_names = narrowPeak_colnames_nostrand)
  
  summit <- (peaks_union$end - peaks_union$start) %/% 2
  peaks_union$summit <- peaks_union$start + summit
  
  #peaks_union %>% write_tsv(., fs::path(path_output_merged, paste0("merged_peaks.hg19_", marker, ".filtered_K27ac_k4me1_k4me3.hq_signal", ".new_summit", ".narrowPeak")),col_names = F)
}


# (Option 2: use point in the middle between multiple summits)


##






