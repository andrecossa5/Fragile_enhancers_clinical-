
library(tidyverse)

SEED <- 4321
set.seed(SEED)
location <- "local" # 'local' or 'hpc'

path_clustered_peaks <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/downstream/deepTools/heatmap.CtIP_and_GRHL_bw.CtIP_and_GRHL_regions.hq_signal.from_mean.clustered.bed")
path_peaks_union_ctip <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/downstream/peaks_union/merged_peaks.hg19_CtIP.filtered_K27ac_k4me1_k4me3.hq_signal.new_summit.narrowPeak")
path_peaks_union_grhl <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/downstream/peaks_union/merged_peaks.hg19_GRHL.filtered_K27ac_k4me1_k4me3.hq_signal.new_summit.narrowPeak")
path_results <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/downstream/")

if(location == "hpc"){
  path_clustered_peaks <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/downstream/deepTools/heatmap.CtIP_and_GRHL_bw.CtIP_and_GRHL_regions.hq_signal.from_mean.clustered.bed")
  path_peaks_union_ctip <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/downstream/peaks_union/merged_peaks.hg19_CtIP.filtered_K27ac_k4me1_k4me3.hq_signal.new_summit.narrowPeak")
  path_peaks_union_grhl <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/downstream/peaks_union/merged_peaks.hg19_GRHL.filtered_K27ac_k4me1_k4me3.hq_signal.new_summit.narrowPeak")
  path_results <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/downstream/")
}

MARKERS <- c("CtIP", "GRHL") 
narrowPeak_colnames_nostrand <- c("chrom", "start", "end", "name", "score", 
                                  "signal_value", "metric1", "metric2", "summit")


##


# Load data
clust_heatmap <- read_tsv(path_clustered_peaks)

peaks_union <- list(
  "CtIP" = read_tsv(path_peaks_union_ctip, col_names = narrowPeak_colnames_nostrand), 
  "GRHL" = read_tsv(path_peaks_union_grhl, col_names = narrowPeak_colnames_nostrand)
)

# All CtIP and GRHL peaks are present in the heatmap 
sum(!peaks_union$CtIP$name %in% clust_heatmap$name)
sum(!peaks_union$GRHL$name %in% clust_heatmap$name)


for(marker in MARKERS){
  
  # Add cluster info based on heatmap
  anno_clust <- peaks_union[[marker]] %>% left_join(., clust_heatmap[,c("name", "deepTools_group")], by = "name")

  print(paste0("Do all ", marker, " enhancers belong to a cluster?"))
  print(sum(is.na(anno_clust$deepTools_group)) == 0)
  
  # Format table for further analyses
  anno_clust <- anno_clust %>% dplyr::select(., chrom, start, end, name, summit, deepTools_group) %>%
    rename(., "deepTools_group" = "cluster")
  
  # Save 
  anno_clust %>% write_tsv(., fs::path(path_results, paste0(marker, "_enh.hq_signal.clustered.tsv")))
  
}


##






