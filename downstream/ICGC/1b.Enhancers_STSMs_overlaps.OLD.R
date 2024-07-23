

library(tidyverse)
library(GenomicRanges)
library(reshape2)
library(reshape2)
library(hrbrthemes) # fro theme_ipsum()
library(moments) # to compute distribution skewness()

SEED <- 4321
set.seed(SEED)

source("/Users/ieo6983/Desktop/fragile_enhancer_clinical/utils/functions_genomics.R")
path_STSMs <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/genomics/pre_processed_ICGC/structural_somatic_mutation.preprocessed.tsv")
path_enhancers_ctip <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/others/Cluster_CtIP_Enh_All.txt")
path_enhancers_grhl <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/others/Cluster_GRHL_Enh_All.txt")


WIN <- 3000
MARKERS <- c("CtIP", "GRHL") 
save_table_overlaps <- F


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


## Load data

STSMs <- read_tsv(path_STSMs)

columns_names <- c("chrom", "start", "end", "cluster")
enh_ctip <- read_tsv(path_enhancers_ctip, col_names = columns_names, comment = "#") 
enh_grhl <- read_tsv(path_enhancers_grhl, col_names = columns_names, comment = "#") 
enh_all <- list("CtIP" = enh_ctip, "GRHL" = enh_grhl)


## Pre-process 

# No start-end, but chrom_bkpt
STSMs <- STSMs[, -c(10, 12, 15, 17)]
STSMs$chrom <- paste("chr", STSMs$chrom, sep = "")
STSMs$start <- STSMs$chrom_bkpt
STSMs$end <- STSMs$chrom_bkpt

# add summit & name to enhancers 
enh_all <- lapply(enh_all, function(df){
  df$summit <- df$end
  df$name <- str_c(df$chrom, df$summit, sep=":")
  return(df)})


## 


# Find overlaps

input_variants <- STSMs
output_enh_STSMs <- list()

for(marker in MARKERS){
  print(marker)
  mut <- "STSMs"
  
  ## Pre/process input data 
  marker_enh <- enh_all[[marker]]
  
  # Extend regions of 'win' from summit
  marker_enh$start <- marker_enh$summit - WIN
  marker_enh$end <- marker_enh$summit + WIN
  
  # Create GRanges objects 
  enh_gr <- makeGRangesFromDataFrame(marker_enh, keep.extra.columns = T)
  STSMs_gr <- makeGRangesFromDataFrame(input_variants, keep.extra.columns = T)
  
  
  ##
  
  
  ## Find mutations / enhancers overlaps 
  hits_obj_enh <- findOverlaps(query=enh_gr, subject=STSMs_gr)
  STSMs_sbj <- data.frame(STSMs_gr[subjectHits(hits_obj_enh)])
  colnames(STSMs_sbj)[1:5] <- paste(colnames(STSMs_sbj)[1:5], "sbj", sep = "_")
  enh_STSMs <- cbind(data.frame(enh_gr[queryHits(hits_obj_enh)]), STSMs_sbj) 
  print( paste0("Total number of overlaps found: ", dim(enh_STSMs)[1], 
                " for window = ", WIN, " bp") )
  
  output_enh_STSMs[[marker]] <- enh_STSMs
  
  
  ##
  
  
  # Save enh_SNVs object for subsequent analyses
  if(save_table_overlaps == T){
    write_tsv(enh_SNVs, 
              fs::path(path_results_data, paste0("Table_enh_SNVs.", marker, ".all_overlaps.", WIN, "bp_WIN", ".tsv")))
  }
}


##


# Add high-low
for(marker in MARKERS){
  output_enh_STSMs[[marker]] <- output_enh_STSMs[[marker]] %>% mutate(., group = ifelse(cluster %in% clust_all$high[[marker]], "high", "low"))
}


for(marker in MARKERS){
  
  # Load data & add high-low 
  table_overlaps <- output_enh_STSMs[[marker]]
  
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
