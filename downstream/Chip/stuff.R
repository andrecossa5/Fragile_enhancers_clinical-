
library(tidyverse)


##


MARKERS <- c("CtIP", "GRHL")
path_input <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/Chip/Chip_for_clusters/results/")
path_input_peaks <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL/macs2/")
path_peaks_union <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL/peaks_union/")
path_output_merged <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL/peaks_union/")


##


# Plots 
num_peaks <- read_tsv("./Desktop/fragile_enhancer_clinical/data/functional_genomics/Chip/Chip_for_clusters/results/macs2/qc/macs2_peak.summary.txt")
num_peaks[num_peaks$measure == "fold", ] %>% 
  ggplot(., aes(x = sample, y = num_peaks, fill = sample))+
  geom_bar(stat = "identity")+
  theme_light()+
  labs(title = "Number of peaks", y = "n. peaks")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "none")
  

anno <- read_tsv("./Desktop/fragile_enhancer_clinical/data/functional_genomics/Chip/Chip_for_clusters/results/macs2/qc/macs2_annotatePeaks.summary.txt")
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

narrowPeak_colnames <- c("chrom", "start", "end", "name", "score", "strand", 
                         "signal_value", "metric1", "metric2", "summit")
peaks <- read_tsv(fs::path(path_input, "CtIP_GRHL/macs2/MCF10A_S30627_NT_201112_hg19_CtIP_Total_peaks.narrowPeak"), 
                  col_names = narrowPeak_colnames)
peaks_anno <- read_tsv(fs::path(path_input, "CtIP_GRHL/macs2/MCF10A_S30627_NT_201112_hg19_CtIP_Total_peaks.annotatePeaks.txt"))
colnames(peaks_anno)[1] <- "PeakID"

# For each .narrowPeak, a txt file with peaks annotation is present
sum(!peaks$name %in% peaks_anno$PeakID)


##


# Filter out potential enhancer peaks 

anno_of_int <- tolower(c("Intergenic", "intron")) # leaving out: TSS, Promoter, exon

peaks_all <- list.files(path_input_peaks)[endsWith(list.files(path_input_peaks), ".narrowPeak")]
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
    
    write_tsv(peaks_bed_filt, fs::path(path_input_peaks, paste0(sample_name, ".filt.narrowPeak" )), col_names = F)
    
  }
}


##
  

# Check if union is correct 
narrowPeak_colnames <- c("chrom", "start", "end", "name", "score", "strand", 
                         "signal_value", "metric1", "metric2", "summit")

ctip_peaks <- read_tsv(fs::path(path_input_peaks, "MCF10A_S30627_NT_201112_hg19_CtIP_Total_peaks.filt.narrowPeak"), col_names = narrowPeak_colnames)
grhl_peaks <- read_tsv(fs::path(path_input_peaks, "MCF10A_S34023_NT_210512_hg19_GRHL2_peaks.filt.narrowPeak"), col_names = narrowPeak_colnames)

ctip_peaks_union <- read_tsv(fs::path(path_peaks_union, "merged_peaks.CtIP.narrowPeak"), col_names = narrowPeak_colnames)
ctip_peaks_union$strand <- NULL
grhl_peaks_union <- read_tsv(fs::path(path_peaks_union, "merged_peaks.GRHL.narrowPeak"), col_names = narrowPeak_colnames)
grhl_peaks_union$strand <- NULL

vec <- findOverlaps(makeGRangesFromDataFrame(ctip_peaks), makeGRangesFromDataFrame(ctip_peaks_union), select = "arbitrary")
vec2 <- findOverlaps(makeGRangesFromDataFrame(grhl_peaks), makeGRangesFromDataFrame(grhl_peaks_union), select = "arbitrary")

sum(is.na(vec))
sum(is.na(vec2))
