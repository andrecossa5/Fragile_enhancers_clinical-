
library(tidyverse)
library(rowr)
library(GenomicRanges)


##


## Section: CHIPs inventory ##
# Comments for which data files to keep or remove, based on quality or different signals
# Specific files for CtIP and GRHL1/2 markers are listed under keep or remove categories.

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


## Section: Paths and settings ##

SEED <- 4321  # Set seed for reproducibility
set.seed(SEED)

# Define markers to be analyzed
MARKERS <- c("CtIP", "GRHL")
location <- "local" # either 'local' or 'hpc' 

# Paths to input and output data files for different ChIP-seq experiments
path_input_peaks <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/macs2/")

# Paths to histone marks (H3K4me3, H3K4me1, H3K27ac)
path_peaks_k4me3 <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/Chip/Chip_for_clusters/results/HISTONE_MARKS/macs2/MCF10A_S48547_230308_hg19_H3K4me3_peaks.narrowPeak")
path_peaks_k4me1 <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/Chip/Chip_for_clusters/results/HISTONE_MARKS/macs2/MCF10A_S47229_SCR_221230_hg19_H3K4me1_peaks.narrowPeak")
path_peaks_k27ac <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/Chip/Chip_for_clusters/results/HISTONE_MARKS/macs2/MCF10A_S43557_SCR_220801_hg19_H3K27ac_peaks.narrowPeak")

# Paths for output results (based on location)
path_output_k4me3 <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/Chip/Chip_for_clusters/results/HISTONE_MARKS/macs2/")

# Update paths if using high-performance computing (hpc) location
if(location == "hpc"){
  path_input_peaks <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/bwa/mergedLibrary/macs2/narrowPeak/")  
  
  path_peaks_k4me3 <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/Chip/Chip_for_clusters/results/HISTONE_MARKS/bwa/mergedLibrary/macs2/narrowPeak/MCF10A_S48547_230308_hg19_H3K4me3_peaks.narrowPeak")
  path_peaks_k4me1 <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/Chip/Chip_for_clusters/results/HISTONE_MARKS/bwa/mergedLibrary/macs2/narrowPeak/MCF10A_S47229_SCR_221230_hg19_H3K4me1_peaks.narrowPeak")
  path_peaks_k27ac <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/Chip/Chip_for_clusters/results/HISTONE_MARKS/bwa/mergedLibrary/macs2/narrowPeak/MCF10A_S43557_SCR_220801_hg19_H3K27ac_peaks.narrowPeak")
  
  path_output_k4me3 <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/Chip/Chip_for_clusters/results/HISTONE_MARKS/bwa/mergedLibrary/macs2/narrowPeak/")
}


# Define column names for narrowPeak files
# `metric1` and `metric2` can represent pValue and qValue respectively
narrowPeak_colnames <- c("chrom", "start", "end", "name", "score", "strand", 
                         "signal_value", "metric1", "metric2", "summit")
narrowPeak_colnames_nostrand <- c("chrom", "start", "end", "name", "score", 
                                  "signal_value", "metric1", "metric2", "summit")


## Section: Plotting peaks summary ##

# Read and plot the total number of peaks for each ChIP sample
num_peaks <- read_tsv(fs::path(path_input_peaks, "qc/macs2_peak.summary.txt")) %>% suppressMessages()
num_peaks[num_peaks$measure == "fold", ] %>% 
  ggplot(., aes(x = sample, y = num_peaks, fill = sample))+
  geom_bar(stat = "identity")+
  theme_light()+
  labs(title = "Number of peaks", y = "n. peaks")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "none")

# Read and plot the annotations of peaks for each ChIP sample
anno <- read_tsv(fs::path(path_input_peaks, "qc/macs2_annotatePeaks.summary.txt")) %>% suppressMessages()
anno[is.na(anno)] <- 0  # Replace NA values with 0
anno <- anno %>% mutate(., tot = Intergenic+TTS+exon+intron+`promoter-TSS`)

# Normalize annotations by total counts and plot the distribution across samples
anno %>% mutate(., across(Intergenic:`promoter-TSS`, ~ . / tot)) %>% 
  dplyr::select(., -tot) %>%
  pivot_longer(-sample) %>% 
  ggplot(., aes(x=sample, y = value, fill = name))+
  geom_bar(stat = "identity", position = "stack")+
  theme_light()+
  labs(title = "Peaks annotation", y = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))


## Section: Filter potential enhancer peaks based on annotations ##

# Define the annotations of interest for potential enhancers
anno_of_int <- tolower(c("Intergenic", "intron")) # Exclude TSS, Promoter, exon

# List all peaks and their annotations in the input directory
peaks_all <- list.files(path_input_peaks)[endsWith(list.files(path_input_peaks), "peaks.narrowPeak")]
peaks_anno_all <- list.files(path_input_peaks)[endsWith(list.files(path_input_peaks), "annotatePeaks.txt")]

# Loop through each marker (CtIP, GRHL) and process the corresponding peaks
for(marker in MARKERS){
  cat("\n")
  print(marker)
  
  # Filter and print annotation files for the current marker
  peaks_anno_marker <- grep(marker, peaks_anno_all, value = T)
  print(peaks_anno_marker)
  
  # Loop through each annotation file and filter peaks based on annotations of interest
  for(anno_file in peaks_anno_marker){
    sample_name <- str_split(anno_file, "\\.")[[1]][1]  # Extract sample name from file name
    
    peaks_anno <- read_tsv(fs::path(path_input_peaks, anno_file)) %>% suppressMessages()
    colnames(peaks_anno)[1] <- "PeakID"  # Rename first column to "PeakID"
    # Simplify annotation column to basic categories
    peaks_anno$simple_annotation <- str_split(peaks_anno$Annotation, pattern = c("-"), simplify = T)[,1] %>% str_split(., pattern=" ", simplify = T) %>% .[,1]
    
    # Select only peaks that could represent enhancers based on the annotations
    peaks_anno_filt <- peaks_anno %>% dplyr::filter(., tolower(simple_annotation) %in% anno_of_int)
    
    # Load corresponding bed file for the sample and filter it based on the annotation filtering
    peaks_bed <- read_tsv(fs::path(path_input_peaks, peaks_all[startsWith(peaks_all, sample_name)]), 
                          col_names = narrowPeak_colnames) %>% suppressMessages() 
    peaks_bed_filt <- peaks_bed[peaks_bed$name %in% peaks_anno_filt$PeakID, ]
    print(dim(peaks_bed_filt)[1] == dim(peaks_anno_filt)[1])  # Check if the number of filtered peaks matches
    
    # Uncomment the next line to save the filtered bed file
    #write_tsv(peaks_bed_filt, fs::path(path_input_peaks, paste0(sample_name, ".filt.narrowPeak" )), col_names = F)
  }
}


##


# --> Merge peaks with bedtools 
# Peaks were subsequently merged using bedtools (command line tool) using script: 3.bedtools.merge_peaks.sh 


##

