
library(tidyverse)
library(GenomicRanges)
library(reshape2)

# Set seed for reproducibility
SEED <- 4321
set.seed(SEED)

location <- "local" # 'local' or 'hpc'

# Load custom genomic functions
source("/Users/ieo6983/Desktop/fragile_enhancer_clinical/utils/functions_genomics.R")

# Define file paths for input data and results
path_SSMs <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/genomics/pre_processed_ICGC/simple_somatic_mutation.open.matching_calls.with_AFs.tsv")
path_enhancers_ctip <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/downstream/CtIP_enh.hq_signal.clustered.tsv")
path_enhancers_grhl <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/downstream/GRHL_enh.hq_signal.clustered.tsv")
path_overlaps_ctip <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/results/ICGC/NEW/enhancers_SSMs_overlaps/data/Table_enh_SNVs.CtIP.all_overlaps.3000bp_WIN.tsv")
path_overlaps_grhl <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/results/ICGC/NEW/enhancers_SSMs_overlaps/data/Table_enh_SNVs.GRHL.all_overlaps.3000bp_WIN.tsv")

path_results_plots <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/results/ICGC/NEW/enhancers_SSMs_overlaps/plots/")

# Override paths and source functions if using HPC environment
if(location == "hpc"){
  source("/hpcnfs/scratch/PGP/Ciacci_et_al/fragile_enhancer_clinical/utils/functions_genomics.R")
  
  path_SSMs <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/data/genomics/pre_processed_ICGC/simple_somatic_mutation.open.matching_calls.with_AFs.tsv")
  path_enhancers_ctip <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/downstream/CtIP_enh.hq_signal.clustered.tsv")
  path_enhancers_grhl <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/downstream/GRHL_enh.hq_signal.clustered.tsv")
  path_overlaps_ctip <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/results/ICGC/NEW/enhancers_SSMs_overlaps/data/Table_enh_SNVs.CtIP.all_overlaps.3000bp_WIN.tsv")
  path_overlaps_grhl <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/results/ICGC/NEW/enhancers_SSMs_overlaps/data/Table_enh_SNVs.GRHL.all_overlaps.3000bp_WIN.tsv")
  
  path_results_plots <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/results/ICGC/NEW/enhancers_SSMs_overlaps/plots/")
}

# Combine file paths into a list for easy access
path_overlaps <- list(
  "CtIP" = path_overlaps_ctip,
  "GRHL" = path_overlaps_grhl
)

WIN <- 3000 # Set window size for enhancer analysis
MARKERS <- c("CtIP", "GRHL") # Set window size for enhancer analysis


##


## High - Low clusters definition

# Clusters are the same for CtIP and GRHL
clust_high_ctip <- c("cluster_1", "cluster_2", "cluster_3", "cluster_4")
clust_low_ctip <- c("cluster_5")

# Create lists of high and low clusters for each marker
clust_high <- list("CtIP" = clust_high_ctip, "GRHL" = clust_high_ctip)
clust_low <- list("CtIP" = clust_low_ctip, "GRHL" = clust_low_ctip)
clust_all <- list("high" = clust_high, "low" = clust_low)


##


# Load overlap data and add high-low group labels based on clusters
table_overlaps_markers <- list(
  "CtIP" = read_tsv(path_overlaps_ctip) %>% 
    mutate(., group = ifelse(cluster %in% clust_all$high$CtIP, "high", "low")), 
  "GRHL" = read_tsv(path_overlaps_grhl) %>% 
    mutate(., group = ifelse(cluster %in% clust_all$high$GRHL, "high", "low"))
) %>% suppressMessages()  # Suppress messages during loading



## Plot distribution of SNVs over enhancers regions

for(marker in MARKERS){
  
  # Extract overlap data for the current marker
  table_overlaps <- table_overlaps_markers[[marker]]
  
  # Calculate the distance of each SNV from the enhancer summit
  dist_enh <- data.frame("dist_enh" = table_overlaps$start_sbj - table_overlaps$summit, 
                         "cluster" = table_overlaps$cluster, 
                         "group" = table_overlaps$group)
  
  ## Plot histogram for all enhancers
  nbins <- 100  # Set number of bins for histogram
  ifelse(marker == "CtIP", color <- "lightblue", color <- "pink")  # Set color based on marker

  ph <- dist_enh %>% ggplot(.)+
    geom_histogram(aes(x = dist_enh), fill = color, alpha = 0.7, bins = nbins)+
    theme_light()+
    scale_x_continuous(breaks = seq(-3000, 3000, 500))+
    labs(title = paste0("Histogram of distances from peak summit - ", marker), 
         subtitle = paste0("Binsize = ", WIN*2/nbins, "bp"),
         x = "distance from summit", y = "")
  print(ph)  
  
  ## Plot density for all enhancers or groups 
  pd <- dist_enh %>% ggplot(.)+
    geom_density(aes(x = dist_enh, fill = group, color = group), alpha = 0.4, bw = nbins)+
    theme_light()+
    scale_x_continuous(breaks = seq(-3000, 3000, 500))+
    labs(title = paste0("Density of distances from peak summit - ", marker), 
         subtitle = paste0("Binsize = ", WIN*2/nbins, "bp"), 
         x = "distance from summit", y = "")
  print(pd) 
  
  # Plot density for each cluster separately
  pd <- dist_enh %>% ggplot(.)+
    geom_density(aes(x = dist_enh, fill = cluster, color = cluster), alpha = 0, bw = nbins)+
    theme_light()+
    scale_x_continuous(breaks = seq(-3000, 3000, 500))+
    labs(title = paste0("Histogram of distances from peak summit - ", marker), 
         subtitle = paste0("Binsize = ", WIN*2/nbins, "bp"),
         x = "distance from summit", y = "")
  print(pd) 
  
  
  # Plot each clusters separately with faceting 
  pc <- dist_enh %>% ggplot(.)+
    geom_density(aes(x = dist_enh, fill = cluster, color = cluster), alpha = 0.7, bw = nbins)+
    theme_light()+
    scale_x_continuous(breaks = seq(-3000, 3000, 500))+
    labs(title = paste0("Histogram of distances from peak summit - ", marker), 
         subtitle = paste0("Binsize = ", WIN*2/nbins, "bp"),
         x = "distance from summit", y = "")+
    facet_wrap(~cluster, ncol = 1)+ # Facet by cluster
    theme(
      strip.text = element_text(size = 5 , margin = margin(1,1,1,1)), # Adjust facet text size
    )
  print(pc) 
  
}


##

