
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
location <- "local" # 'local' or 'hpc'

# Paths CtIP/GRHL peaks
path_input_peaks <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/macs2/")

# Paths histone marks 
path_peaks_k4me3 <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/Chip/Chip_for_clusters/results/HISTONE_MARKS/macs2/MCF10A_S48547_230308_hg19_H3K4me3_peaks.narrowPeak")
path_peaks_k4me1 <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/Chip/Chip_for_clusters/results/HISTONE_MARKS/macs2/MCF10A_S47229_SCR_221230_hg19_H3K4me1_peaks.narrowPeak")
path_peaks_k27ac <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/Chip/Chip_for_clusters/results/HISTONE_MARKS/macs2/MCF10A_S43557_SCR_220801_hg19_H3K27ac_peaks.narrowPeak")

# Paths output 
path_output_k4me3 <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/Chip/Chip_for_clusters/results/HISTONE_MARKS/macs2/")

if(location == "hpc"){
  path_input_peaks <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/bwa/mergedLibrary/macs2/narrowPeak/")  
  
  path_peaks_k4me3 <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/Chip/Chip_for_clusters/results/HISTONE_MARKS/bwa/mergedLibrary/macs2/narrowPeak/MCF10A_S48547_230308_hg19_H3K4me3_peaks.narrowPeak")
  path_peaks_k4me1 <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/Chip/Chip_for_clusters/results/HISTONE_MARKS/bwa/mergedLibrary/macs2/narrowPeak/MCF10A_S47229_SCR_221230_hg19_H3K4me1_peaks.narrowPeak")
  path_peaks_k27ac <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/Chip/Chip_for_clusters/results/HISTONE_MARKS/bwa/mergedLibrary/macs2/narrowPeak/MCF10A_S43557_SCR_220801_hg19_H3K27ac_peaks.narrowPeak")
  
  path_output_k4me3 <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/Chip/Chip_for_clusters/results/HISTONE_MARKS/bwa/mergedLibrary/macs2/narrowPeak/")
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


##

