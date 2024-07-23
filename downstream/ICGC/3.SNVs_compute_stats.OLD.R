
library(tidyverse)
library(GenomicRanges)
library(reshape2)
library(rowr)
library(ggpubr)

SEED <- 4321
set.seed(SEED)

source("/Users/ieo6983/Desktop/fragile_enhancer_clinical/utils/functions_genomics.R")

path_SSMs <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/genomics/pre_processed_ICGC/simple_somatic_mutation.open.matching_calls.with_AFs.tsv")
path_enhancers_ctip <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/others/Cluster_CtIP_Enh_All.txt")
path_enhancers_grhl <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/others/Cluster_GRHL_Enh_All.txt")

path_overlaps_ctip <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/results/ICGC/OLD/enhancers_SSMs_overlaps/data/Table_enh_SNVs.CtIP.all_overlaps.3000bp_WIN.tsv")
path_overlaps_grhl <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/results/ICGC/OLD/enhancers_SSMs_overlaps/data/Table_enh_SNVs.GRHL.all_overlaps.3000bp_WIN.tsv")

path_results_data <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/results/ICGC/OLD/enhancers_SSMs_overlaps/data/")
path_results_plots <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/results/ICGC/OLD/enhancers_SSMs_overlaps/plots/")

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


## Load input files 
SSMs <- read_tsv(path_SSMs)

columns_names <- c("chrom", "start", "end", "cluster")
enh_ctip <- read_tsv(path_enhancers_ctip, col_names = columns_names, comment = "#") 
enh_grhl <- read_tsv(path_enhancers_grhl, col_names = columns_names, comment = "#") 
enh_all <- list("CtIP" = enh_ctip, "GRHL" = enh_grhl)


## Pre-process 

# add summit & name to enhancers 
enh_all <- lapply(enh_all, function(df){
  df$summit <- df$end
  df$name <- str_c(df$chrom, df$summit, sep=":")
  return(df)})
for(marker in MARKERS){
  enh_all[[marker]]$group <- ifelse(enh_all[[marker]]$cluster %in% clust_all$high[[marker]], "high", "low")
}

# Filter SSMs (only SNVs)
SNVs <- SSMs %>% rename(., c("chromosome_start" = "start", "chromosome_end" = "end")) %>% 
  mutate(., chromosome = paste0("chr", .$chromosome)) %>%
  dplyr::filter(., end - start == 0)

# Load overlaps & add high-low labels
table_overlaps_markers <- list(
  "CtIP" = read_tsv(path_overlaps_ctip) %>% mutate(., group = ifelse(cluster %in% clust_high_ctip, "high", "low")), 
  "GRHL" = read_tsv(path_overlaps_grhl) %>% mutate(., group = ifelse(cluster %in% clust_high_grhl, "high", "low"))
) %>% suppressMessages()


##


# Find overlaps with random regions 

input_variants <- SNVs
ran_seqs_df <- list()
ran_seqs_overlaps <- list()

for(marker in MARKERS){
  print(marker)
  marker_enh <- enh_all[[marker]]
  
  SNVs_gr <- makeGRangesFromDataFrame(input_variants, keep.extra.columns = T)
  
  # Generate random sequences 
  ran_seqs <- generate_random_seqs(SEED, dim(marker_enh)[1], win = WIN)
  ran_seqs$summit <- ran_seqs$start + (WIN/2)
  ran_seqs_df[[marker]] <- ran_seqs
  
  # Compute overlaps 
  ran_seqs_gr = makeGRangesFromDataFrame(ran_seqs, keep.extra.columns = T)
  
  hits_obj_ran <- findOverlaps(query=ran_seqs_gr, subject=SNVs_gr)
  SNVs_sbj <- data.frame(SNVs_gr[subjectHits(hits_obj_ran)])
  colnames(SNVs_sbj)[1:5] <- paste(colnames(SNVs_sbj)[1:5], "sbj", sep = "_")
  ran_seqs_SNVs <- cbind(data.frame(ran_seqs_gr[queryHits(hits_obj_ran)], SNVs_sbj))
  
  ran_seqs_overlaps[[marker]] <- ran_seqs_SNVs
  
}


##


## Compute stats for multiple windows around summits

WINS <- c(3000, 2000, 1000, 500, 100, 50)

stat_enh_all_wins <- list()
stat_vars_all_wins <- list()

suppressMessages({

for(marker in MARKERS){
  
  for(win in WINS){
    
    overlaps <- table_overlaps_markers[[marker]]
    input_enh <- enh_all[[marker]]
    
    # Extract overlaps respecting window
    overlaps$dist <- abs(overlaps$start_sbj - overlaps$summit)
    overlaps <- overlaps %>% dplyr::filter(., dist <= win)
  
    
    ## Compute stats for enhancers
    
    # Number of donors in which each enhancer is hit
    stat_enh <- overlaps %>% dplyr::group_by(name, group) %>%
      dplyr::summarise(., n_donors_hit = length(unique(icgc_donor_id)))
    # Add to df enhancers with no mutation in any patient
    stat_enh <- rbind(stat_enh, data.frame("name" = input_enh$name[!input_enh$name %in% stat_enh$name], 
                                           "group" = input_enh[!input_enh$name %in% stat_enh$name, ]$group,
                                           "n_donors_hit" = rep(0, length(input_enh$name[!input_enh$name %in% stat_enh$name]))))
    
    # Number of variants within enhancer 
    stat_vars <- overlaps %>% group_by(name, group) %>% 
      dplyr::summarise(., n_vars = n())
    stat_vars <- rbind(stat_vars, data.frame("name" = input_enh$name[!input_enh$name %in% stat_vars$name], 
                                             "group" = input_enh[!input_enh$name %in% stat_vars$name, ]$group,
                                           "n_vars" = rep(0, length(input_enh$name[!input_enh$name %in% stat_vars$name]))))
    
    
    ##
    
    
    ## Compute stats for random seqs
    overlaps <- ran_seqs_overlaps[[marker]]
    input_enh <- ran_seqs_df[[marker]]
    
    # Extract overlaps respecting window
    overlaps$dist <- abs(overlaps$start_sbj - overlaps$summit)
    overlaps <- overlaps %>% dplyr::filter(., dist <= win)
    
    # Number of donors in which each enhancer is hit
    stat_enh_ran <- overlaps %>% dplyr::group_by(name) %>%
      dplyr::summarise(., n_donors_hit = length(unique(icgc_donor_id)))
    # Add to df enhancers with no mutation in any patient
    stat_enh_ran <- rbind(stat_enh_ran, data.frame("name" = input_enh$name[!input_enh$name %in% stat_enh_ran$name], 
                                               "n_donors_hit" = rep(0, length(input_enh$name[!input_enh$name %in% stat_enh_ran$name]))))
    # Number of variants within enhancer 
    stat_vars_ran <- overlaps %>% group_by(name) %>% 
      dplyr::summarise(., n_vars = n())
    stat_vars_ran <- rbind(stat_vars_ran, data.frame("name" = input_enh$name[!input_enh$name %in% stat_vars_ran$name], 
                                             "n_vars" = rep(0, length(input_enh$name[!input_enh$name %in% stat_vars_ran$name]))))
    
    # Merge & save 
    stat_enh_all <- cbind.fill(stat_enh, stat_enh_ran)
    colnames(stat_enh_all) <- c("name", "group", "n_donors_hit", "name.ran", "n_donors_hit_ran")
    stat_enh_all$n_donors_hit <- as.numeric(stat_enh_all$n_donors_hit)
    stat_enh_all$n_donors_hit_ran <- as.numeric(stat_enh_all$n_donors_hit_ran)
    
    stat_vars_all <- cbind.fill(stat_vars, stat_vars_ran)
    colnames(stat_vars_all) <- c("name", "group", "n_vars", "name.ran", "n_vars_ran")
    stat_vars_all$n_vars <- as.numeric(stat_vars_all$n_vars)
    stat_vars_all$n_vars_ran <- as.numeric(stat_vars_all$n_vars_ran)
    
    win_c <- as.character(win)
    stat_enh_all_wins[[marker]][[win_c]] <- stat_enh_all
    stat_vars_all_wins[[marker]][[win_c]] <- stat_vars_all
    
  }
}

})


##
  

## Plot stats
for(marker in MARKERS){
  for(win in WINS){
    win_c <- as.character(win)
    stat_enh_to_plot <- stat_enh_all_wins[[marker]][[win_c]]
    stat_vars_to_plot <- stat_vars_all_wins[[marker]][[win_c]]
    
    # Compute mean for each group 
    avg <- stat_enh_to_plot %>% group_by(., group) %>% dplyr::summarise(., avg = mean(n_donors_hit))
    avg <- str_flatten(paste0(avg$group, ": ", round(avg$avg,3)), collapse = " ")
  
    p1 <- stat_enh_to_plot %>% dplyr::select(., n_donors_hit, n_donors_hit_ran, group) %>%
      pivot_longer(-group) %>%
      ggplot(., aes(x = group, y = value, fill = name))+
      geom_violin()+
      labs(title = paste0("Number of donors hit - WIN ", win), 
           subtitle = avg,
           x = "", y = "number of donors")+
      theme_light()+
      stat_compare_means(method = "wilcox", label = "p.signif", 
                         comparisons = list( c("high", "low")))+
      stat_compare_means(method = "wilcox", label = "p.signif")
    print(p1)
    
    
    avg <- stat_vars_to_plot %>% group_by(., group) %>% dplyr::summarise(., avg = mean(n_vars))
    avg <- str_flatten(paste0(avg$group, ": ", round(avg$avg,3)), collapse = " ")
    
    p2 <- stat_vars_to_plot %>% dplyr::select(., n_vars, n_vars_ran, group) %>%
      pivot_longer(-group) %>%
      ggplot(., aes(x = group, y = value, fill = name))+
      geom_violin()+
      labs(title = paste0("Number of variants x enhancer - WIN ", win),
           subtitle = avg,
           x = "", y = "number of donors")+
      theme_light()+
      stat_compare_means(method = "wilcox", label = "p.signif", 
                         comparisons = list( c("high", "low")))+
      stat_compare_means(method = "wilcox", label = "p.signif")
    print(p2)
  }
}


##


