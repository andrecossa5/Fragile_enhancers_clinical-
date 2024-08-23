
suppressMessages({
  
library(tidyverse)
library(GenomicRanges)
library(reshape2)
library(rowr)
library(ggpubr)

# Set seed for reproducibility
SEED <- 4321
set.seed(SEED)

# Source additional functions
source("/hpcnfs/scratch/P_PGP_FRAGILE_ENHANCERS/bin/functions_genomics.R")

# Define file paths
path_SSMs <- fs::path("/hpcnfs/scratch/P_PGP_FRAGILE_ENHANCERS/results/data/Hartwig_all_snvs_info.tsv")
path_enhancers_ctip <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/downstream/CtIP_enh.hq_signal.clustered.tsv")
path_enhancers_grhl <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/downstream/GRHL_enh.hq_signal.clustered.tsv")
path_overlaps_ctip <- fs::path("/hpcnfs/scratch/P_PGP_FRAGILE_ENHANCERS/results/NEW/data/CtIP_enh.hartwig_snvs.overlap.WIN_3000.tsv")
path_overlaps_grhl <- fs::path("/hpcnfs/scratch/P_PGP_FRAGILE_ENHANCERS/results/NEW/data/GRHL_enh.hartwig_snvs.overlap.WIN_3000.tsv")
path_results_plots <- fs::path("/hpcnfs/scratch/P_PGP_FRAGILE_ENHANCERS/results/NEW/plots")

# Define parameters
WIN <- 3000
MARKERS <- c("CtIP", "GRHL")


##

# Define the file name and open PDF device for plotting
file_name <- fs::path(path_results_plots, paste0("Rplots.stats.snvs_frequency_enhancers_and_random_controls.pdf")) 
pdf(file_name)

##


## High - Low clusters definition

# Define high and low clusters for CtIP and GRHL
clust_high_ctip <- c("cluster_1", "cluster_2", "cluster_3", "cluster_4")
clust_low_ctip <- c("cluster_5")

# Assign high and low clusters to markers
clust_high <- list("CtIP" = clust_high_ctip, "GRHL" = clust_high_ctip)
clust_low <- list("CtIP" = clust_low_ctip, "GRHL" = clust_low_ctip)
clust_all <- list("high" = clust_high, "low" = clust_low)


##


## Load input files 

# Read data for SNVs and enhancers
SSMs <- read_tsv(path_SSMs)
enh_ctip <- read_tsv(path_enhancers_ctip) 
enh_grhl <- read_tsv(path_enhancers_grhl) 
enh_all <- list("CtIP" = enh_ctip, "GRHL" = enh_grhl)


## Pre-process 

# Modify enhancer names and group them into high and low clusters
enh_all <- lapply(enh_all, function(df){
  df$chrom <- paste0("chr", df$chrom)
  df$name <- str_c(df$chrom, df$summit, sep=":")
  return(df)})

# Assign high or low group based on clusters
for(marker in MARKERS){
  enh_all[[marker]]$group <- ifelse(enh_all[[marker]]$cluster %in% clust_all$high[[marker]], "high", "low")
}

# Filter SNVs to include only single nucleotide variants
SNVs <- SSMs %>% rename(., c("seqnames" = "chromosome")) %>% 
  dplyr::filter(., width == 1)

# Load overlaps data and add high-low labels
column_names <- c(
  paste0(c("seqnames", "start", "end", "width", "strand", "name", "summit", "cluster"), "_enh"), 
  paste0(c("seqnames", "start", "end", "width", "strand", "REF", "ALT", "QUAL", "FILTER", "SAMPLE", "PURPLE_AF", "AF", "ID"), "_snv")
)
table_overlaps_markers <- list(
  "CtIP" = read_tsv(path_overlaps_ctip), 
  "GRHL" = read_tsv(path_overlaps_grhl)
)
for(marker in MARKERS){
  colnames(table_overlaps_markers[[marker]]) <- column_names  # Rename columns
  table_overlaps_markers[[marker]] <- table_overlaps_markers[[marker]] %>% mutate(., group = ifelse(cluster_enh %in% clust_all$high[[marker]], "high", "low"))
}


##


# Find overlaps with random regions 

input_variants <- SNVs
ran_seqs_df <- list()
ran_seqs_overlaps <- list()

for(marker in MARKERS){
  print(marker)
  marker_enh <- enh_all[[marker]]
  
  # Convert SNVs to GRanges object
  SNVs_gr <- makeGRangesFromDataFrame(input_variants, keep.extra.columns = T)
  
  # Generate random sequences 
  ran_seqs <- generate_random_seqs(SEED, dim(marker_enh)[1], win = WIN)
  ran_seqs$summit <- ran_seqs$start + (WIN/2)
  ran_seqs_df[[marker]] <- ran_seqs
  
  # Compute overlaps between random sequences and SNVs
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

for(marker in MARKERS){
  
  for(win in WINS){
    
    overlaps <- table_overlaps_markers[[marker]]
    input_enh <- enh_all[[marker]]
    
    # Extract overlaps within the window size
    overlaps$dist <- abs(overlaps$start_snv - overlaps$summit_enh)
    overlaps <- overlaps %>% dplyr::filter(., dist <= win)
  
    
    ## Compute stats for enhancers
    
    # Number of donors in which each enhancer is hit
    stat_enh <- overlaps %>% dplyr::group_by(name_enh, group) %>%
      dplyr::summarise(., n_donors_hit = length(unique(SAMPLE_snv)))
    # Add to df enhancers with no mutation in any patient
    stat_enh <- rbind(stat_enh, data.frame("name_enh" = input_enh$name[!input_enh$name %in% stat_enh$name_enh], 
                                           "group" = input_enh[!input_enh$name %in% stat_enh$name_enh, ]$group,
                                           "n_donors_hit" = rep(0, length(input_enh$name[!input_enh$name %in% stat_enh$name_enh]))))
    
    # Number of variants within enhancer 
    stat_vars <- overlaps %>% group_by(name_enh, group) %>% 
      dplyr::summarise(., n_vars = n())
    stat_vars <- rbind(stat_vars, data.frame("name_enh" = input_enh$name[!input_enh$name %in% stat_vars$name_enh], 
                                             "group" = input_enh[!input_enh$name %in% stat_vars$name_enh, ]$group,
                                           "n_vars" = rep(0, length(input_enh$name[!input_enh$name %in% stat_vars$name_enh]))))
    
    
    ##
    
    
    ## Compute stats for random sequences
    overlaps <- ran_seqs_overlaps[[marker]]
    input_enh <- ran_seqs_df[[marker]]
    
    # Extract overlaps within the window size
    overlaps$dist <- abs(overlaps$start_sbj - overlaps$summit)
    overlaps <- overlaps %>% dplyr::filter(., dist <= win)
    
    # Count the number of donors hitting each enhancer in random sequences
    stat_enh_ran <- overlaps %>% dplyr::group_by(name) %>%
      dplyr::summarise(., n_donors_hit = length(unique(SAMPLE)))
    # Add enhancers with no mutations in any patient
    stat_enh_ran <- rbind(stat_enh_ran, data.frame("name" = input_enh$name[!input_enh$name %in% stat_enh_ran$name], 
                                                   "n_donors_hit" = rep(0, length(input_enh$name[!input_enh$name %in% stat_enh_ran$name]))))
    # Count the number of variants within each enhancer in random sequences
    stat_vars_ran <- overlaps %>% group_by(name) %>% 
      dplyr::summarise(., n_vars = n())
    stat_vars_ran <- rbind(stat_vars_ran, data.frame("name" = input_enh$name[!input_enh$name %in% stat_vars_ran$name], 
                                                     "n_vars" = rep(0, length(input_enh$name[!input_enh$name %in% stat_vars_ran$name]))))
    
    # Merge and save stats
    stat_enh_all <- cbind.fill(stat_enh, stat_enh_ran)
    colnames(stat_enh_all) <- c("name", "group", "n_donors_hit", "name_ran", "n_donors_hit_ran")
    stat_enh_all$n_donors_hit <- as.numeric(stat_enh_all$n_donors_hit)
    stat_enh_all$n_donors_hit_ran <- as.numeric(stat_enh_all$n_donors_hit_ran)
    
    stat_vars_all <- cbind.fill(stat_vars, stat_vars_ran)
    colnames(stat_vars_all) <- c("name", "group", "n_vars", "name_ran", "n_vars_ran")
    stat_vars_all$n_vars <- as.numeric(stat_vars_all$n_vars)
    stat_vars_all$n_vars_ran <- as.numeric(stat_vars_all$n_vars_ran)
    
    win_c <- as.character(win)
    stat_enh_all_wins[[marker]][[win_c]] <- stat_enh_all
    stat_vars_all_wins[[marker]][[win_c]] <- stat_vars_all
    
  }
}


##
  

## Plot stats

# Create plots for each marker and window size
for(marker in MARKERS){
  for(win in WINS){
    win_c <- as.character(win)
    stat_enh_to_plot <- stat_enh_all_wins[[marker]][[win_c]]
    stat_vars_to_plot <- stat_vars_all_wins[[marker]][[win_c]]
    
    # Compute mean number of donors hit for each group
    avg <- stat_enh_to_plot %>% group_by(., group) %>% dplyr::summarise(., avg = mean(n_donors_hit))
    avg <- str_flatten(paste0(avg$group, ": ", round(avg$avg,3)), collapse = " ")
  
    mean_all <- stat_enh_to_plot %>% dplyr::select(., n_donors_hit, n_donors_hit_ran, group) %>%
      pivot_longer(-group) %>% group_by(., group, name) %>%
      dplyr::summarise(., mean = mean(value))
    
    # Plot number of donors hit
    p1 <- stat_enh_to_plot %>% dplyr::select(., n_donors_hit, n_donors_hit_ran, group) %>%
      pivot_longer(-group) %>%
      ggplot(., aes(x = group, y = value, fill = name))+
      geom_violin()+
      geom_point(data = mean_all, aes(x = group, y = mean, color = name), 
                 shape = 23, size = 3, fill = "white") +
      labs(title = paste0("Number of donors hit - WIN ", win, " - Marker: ", marker), 
           subtitle = avg,
           x = "", y = "number of donors")+
      theme_light()+
      stat_compare_means(method = "wilcox", label = "p.signif", 
                         comparisons = list( c("high", "low")))+
      stat_compare_means(method = "wilcox", label = "p.signif")
    print(p1)
    
    
    # Compute mean number of variants for each group
    avg <- stat_vars_to_plot %>% group_by(., group) %>% dplyr::summarise(., avg = mean(n_vars))
    avg <- str_flatten(paste0(avg$group, ": ", round(avg$avg,3)), collapse = " ")
    
    mean_all <- stat_vars_to_plot %>% dplyr::select(., n_vars, n_vars_ran, group) %>%
      pivot_longer(-group) %>% group_by(., group, name) %>%
      dplyr::summarise(., mean = mean(value))
    
    # Plot number of variants per enhancer
    p2 <- stat_vars_to_plot %>% dplyr::select(., n_vars, n_vars_ran, group) %>%
      pivot_longer(-group) %>%
      ggplot(., aes(x = group, y = value, fill = name))+
      geom_violin()+
      geom_point(data = mean_all, aes(x = group, y = mean, color = name), 
                 shape = 23, size = 3, fill = "white") +
      labs(title = paste0("Number of variants x enhancer - WIN ", win, " - Marker: ", marker),
           subtitle = avg,
           x = "", y = "number of donors")+
      theme_light()+
      stat_compare_means(method = "wilcox", label = "p.signif", 
                         comparisons = list(c("high", "low")))+
      stat_compare_means(method = "wilcox", label = "p.signif")
    print(p2)
  }
}


##

dev.off()

##

})
