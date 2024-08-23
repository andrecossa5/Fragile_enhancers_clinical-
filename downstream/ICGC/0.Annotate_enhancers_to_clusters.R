
library(tidyverse)

# Set the seed for reproducibility
SEED <- 4321
set.seed(SEED)

location <- "local" # 'local' or 'hpc'

# Define file paths for local environment
path_clustered_peaks <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/downstream/deepTools/heatmap.CtIP_and_GRHL_bw.CtIP_and_GRHL_regions.hq_signal.from_mean.clustered.bed")
path_peaks_union_ctip <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/downstream/peaks_union/merged_peaks.hg19_CtIP.filtered_K27ac_k4me1_k4me3.hq_signal.new_summit.narrowPeak")
path_peaks_union_grhl <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/downstream/peaks_union/merged_peaks.hg19_GRHL.filtered_K27ac_k4me1_k4me3.hq_signal.new_summit.narrowPeak")
path_results <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/downstream/")

# Override file paths if using HPC environment
if(location == "hpc"){
  path_clustered_peaks <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/downstream/deepTools/heatmap.CtIP_and_GRHL_bw.CtIP_and_GRHL_regions.hq_signal.from_mean.clustered.bed")
  path_peaks_union_ctip <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/downstream/peaks_union/merged_peaks.hg19_CtIP.filtered_K27ac_k4me1_k4me3.hq_signal.new_summit.narrowPeak")
  path_peaks_union_grhl <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/downstream/peaks_union/merged_peaks.hg19_GRHL.filtered_K27ac_k4me1_k4me3.hq_signal.new_summit.narrowPeak")
  path_results <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/downstream/")
}

# Define the markers to analyze
MARKERS <- c("CtIP", "GRHL") 
# Define column names for narrowPeak files without strand information
narrowPeak_colnames_nostrand <- c("chrom", "start", "end", "name", "score", 
                                  "signal_value", "metric1", "metric2", "summit")


##


## Data Loading ##

# Load the clustered peaks data (heatmap) into a dataframe
clust_heatmap <- read_tsv(path_clustered_peaks)

# Load the peaks data for both CtIP and GRHL into a list of dataframes
peaks_union <- list(
  "CtIP" = read_tsv(path_peaks_union_ctip, col_names = narrowPeak_colnames_nostrand), 
  "GRHL" = read_tsv(path_peaks_union_grhl, col_names = narrowPeak_colnames_nostrand)
)

# Ensure that all CtIP and GRHL peaks are present in the heatmap data
sum(!peaks_union$CtIP$name %in% clust_heatmap$name)
sum(!peaks_union$GRHL$name %in% clust_heatmap$name)

# Loop through each marker (CtIP and GRHL)
for(marker in MARKERS){
  
  # Add cluster information to each peak by joining with the heatmap data
  anno_clust <- peaks_union[[marker]] %>% 
    left_join(., clust_heatmap[,c("name", "deepTools_group")], by = "name")
  
  # Check if all peaks for the current marker are assigned to a cluster
  print(paste0("Do all ", marker, " enhancers belong to a cluster?"))
  print(sum(is.na(anno_clust$deepTools_group)) == 0)
  
  # Select and rename relevant columns for further analysis
  anno_clust <- anno_clust %>% 
    dplyr::select(., chrom, start, end, name, summit, deepTools_group) %>%
    rename(., "deepTools_group" = "cluster")
  
  # Save the annotated peaks to a file
  anno_clust %>% write_tsv(., fs::path(path_results, paste0(marker, "_enh.hq_signal.clustered.tsv")))
  
}


##

