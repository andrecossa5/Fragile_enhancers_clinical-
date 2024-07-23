
library(tidyverse)
library(GenomicRanges)
library(reshape2)
library(reshape2)
library(hrbrthemes) # fro theme_ipsum()
library(moments) # to compute distribution skewness()

SEED <- 4321
set.seed(SEED)

source("/Users/ieo6983/Desktop/fragile_enhancer_clinical/utils/functions_genomics.R")

path_SSMs <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/genomics/pre_processed_ICGC/simple_somatic_mutation.open.matching_calls.with_AFs.tsv")
path_enhancers_ctip <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/others/Cluster_CtIP_Enh_All.txt")
path_enhancers_grhl <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/others/Cluster_GRHL_Enh_All.txt")

path_overlaps_ctip <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/results/ICGC/OLD/enhancers_SSMs_overlaps/data/Table_enh_SNVs.CtIP.all_overlaps.3000bp_WIN.tsv")
path_overlaps_grhl <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/results/ICGC/OLD/enhancers_SSMs_overlaps/data/Table_enh_SNVs.GRHL.all_overlaps.3000bp_WIN.tsv")

path_overlaps <- list(
  "CtIP" = path_overlaps_ctip,
  "GRHL" = path_overlaps_grhl
)

WIN <- 3000
MARKERS <- c("CtIP", "GRHL") 


##


## High - Low clusters definition

# Option 1: High = first 3 clusters; Low = last 4 clusters
clust_high_ctip <- c("CtIP_cluster_2_Enh", "CtIP_cluster_3_Enh", "CtIP_cluster_5_Enh")
clust_high_grhl <- c("GRHL_cluster_1_Enh", "GRHL_cluster_2_Enh", "GRHL_cluster_4_Enh")
clust_low_ctip <- c("CtIP_cluster_6.2_Enh", "CtIP_cluster_6.3_Enh", "CtIP_cluster_6.1_Enh", "CtIP_cluster_6.0")
clust_low_grhl <- c("GRHL_cluster_5.2_Enh", "GRHL_cluster_5.3_Enh", "GRHL_cluster_5.1_Enh", "GRHL_cluster_5.0")

clust_high <- list("CtIP" = clust_high_ctip, "GRHL" = clust_high_grhl)
clust_low <- list("CtIP" = clust_low_ctip, "GRHL" = clust_low_grhl)
clust_all <- list("high" = clust_high, "low" = clust_low)


##


# Load data & add high-low labels
table_overlaps_markers <- list(
  "CtIP" = read_tsv(path_overlaps_ctip) %>% mutate(., group = ifelse(cluster %in% clust_high_ctip, "high", "low")), 
  "GRHL" = read_tsv(path_overlaps_grhl) %>% mutate(., group = ifelse(cluster %in% clust_high_grhl, "high", "low"))
) %>% suppressMessages()


## Plot distribution of SNVs over enhancers regions

for(marker in MARKERS){
  
  # Load data & add high-low 
  table_overlaps <- table_overlaps_markers[[marker]]
  
  # Compute dist of SNVs from enhancer summit
  dist_enh <- data.frame("dist_enh" = table_overlaps$start_sbj - table_overlaps$summit, 
                         "cluster" = table_overlaps$cluster, 
                         "group" = table_overlaps$group)
  
  
  ## Plot histogram for all enhancers
  nbins <- 100
  ifelse(marker == "CtIP", color <- "lightblue", color <- "pink")

  ph <- dist_enh %>% ggplot(.)+
    geom_histogram(aes(x = dist_enh), fill = color, alpha = 0.7, bins = nbins)+
    theme_light()+
    scale_x_continuous(breaks = seq(-3000, 3000, 500))+
    labs(title = paste("Histogram of distances from peak summit - ", marker, " (Binsize = ", WIN*2/nbins, "bp)", sep =""), 
         x = "distance from summit", y = "")
  print(ph)  
  
  ## Plot density for all enhancers or groups 
  pd <- dist_enh %>% ggplot(.)+
    geom_density(aes(x = dist_enh, fill = group, color = group), alpha = 0.4, bw = nbins)+
    theme_light()+
    scale_x_continuous(breaks = seq(-3000, 3000, 500))+
    labs(title = paste("Density of distances from peak summit - ", marker, " (Binsize = ", WIN*2/nbins, "bp)", sep =""), 
         x = "distance from summit", y = "")
  print(pd) 
  
  pd <- dist_enh %>% ggplot(.)+
    geom_density(aes(x = dist_enh, fill = cluster, color = cluster), alpha = 0, bw = nbins)+
    theme_light()+
    scale_x_continuous(breaks = seq(-3000, 3000, 500))+
    labs(title = paste("Density of distances from peak summit - ", marker, " (Binsize = ", WIN*2/nbins, "bp)", sep =""), 
         x = "distance from summit", y = "")
  print(pd) 
  
  
  # Plot each clusters separately 
  pc <- dist_enh %>% ggplot(.)+
    geom_density(aes(x = dist_enh, fill = cluster, color = cluster), alpha = 0.7, bw = nbins)+
    theme_light()+
    scale_x_continuous(breaks = seq(-3000, 3000, 500))+
    labs(title = paste("Density of distances from peak summit - ", marker, " (Binsize = ", WIN*2/nbins, "bp)", sep =""), 
         x = "distance from summit", y = "")+
    facet_wrap(~cluster, ncol = 1)+
    theme(
      strip.text = element_text(size = 5 , margin = margin(1,1,1,1)),
    )
  print(pc) 
}


##





