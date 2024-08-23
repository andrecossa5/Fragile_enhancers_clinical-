
library(tidyverse)
library(GenomicRanges)
library(reshape2)
library(rowr)
library(ggpubr)

# Set random seed for reproducibility
SEED <- 4321
set.seed(SEED)

location <- "local" # 'local' or 'hpc'

# Load custom functions from external file
source("/Users/ieo6983/Desktop/fragile_enhancer_clinical/utils/functions_genomics.R")

# Define file paths for local environment
path_SSMs <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/genomics/pre_processed_ICGC/simple_somatic_mutation.open.matching_calls.with_AFs.tsv")
path_enhancers_ctip <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/downstream/CtIP_enh.hq_signal.clustered.tsv")
path_enhancers_grhl <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/downstream/GRHL_enh.hq_signal.clustered.tsv")
path_overlaps_ctip <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/results/ICGC/NEW/enhancers_SSMs_overlaps/data/Table_enh_SNVs.CtIP.all_overlaps.3000bp_WIN.tsv")
path_overlaps_grhl <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/results/ICGC/NEW/enhancers_SSMs_overlaps/data/Table_enh_SNVs.GRHL.all_overlaps.3000bp_WIN.tsv")

path_results_plots <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/results/ICGC/NEW/enhancers_SSMs_overlaps/plots/")


# If running on HPC, update paths accordingly
if(location == "hpc"){
  source("/hpcnfs/scratch/PGP/Ciacci_et_al/fragile_enhancer_clinical/utils/functions_genomics.R")
  path_SSMs <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/data/genomics/pre_processed_ICGC/simple_somatic_mutation.open.matching_calls.with_AFs.tsv")
  path_enhancers_ctip <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/downstream/CtIP_enh.hq_signal.clustered.tsv")
  path_enhancers_grhl <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/downstream/GRHL_enh.hq_signal.clustered.tsv") 
  path_overlaps_ctip <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/results/ICGC/NEW/enhancers_SSMs_overlaps/data/Table_enh_SNVs.CtIP.all_overlaps.3000bp_WIN.tsv")
  path_overlaps_grhl <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/results/ICGC/NEW/enhancers_SSMs_overlaps/data/Table_enh_SNVs.GRHL.all_overlaps.3000bp_WIN.tsv")
  
  path_results_plots <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/results/ICGC/NEW/enhancers_SSMs_overlaps/plots/")
}

# Define constants
WIN <- 3000
MARKERS <- c("CtIP", "GRHL") 


##


## High - Low clusters definition

# Define high and low clusters for CtIP and GRHL markers - they are the same for CtIP and GRHL
clust_high_ctip <- c("cluster_1", "cluster_2", "cluster_3", "cluster_4")
clust_low_ctip <- c("cluster_5")

# Create lists to map clusters to markers
clust_high <- list("CtIP" = clust_high_ctip, "GRHL" = clust_high_ctip)
clust_low <- list("CtIP" = clust_low_ctip, "GRHL" = clust_low_ctip)
clust_all <- list("high" = clust_high, "low" = clust_low)


##


## Load input files 

# Load somatic mutation data
SSMs <- read_tsv(path_SSMs)

# Load enhancer data for CtIP and GRHL markers
enh_ctip <- read_tsv(path_enhancers_ctip) 
enh_grhl <- read_tsv(path_enhancers_grhl) 
enh_all <- list("CtIP" = enh_ctip, "GRHL" = enh_grhl)


## Pre-process 

# Pre-process enhancer data: add 'chr' prefix to chromosome names and create unique names for each enhancer
enh_all <- lapply(enh_all, function(df){
  df$chrom <- paste0("chr", df$chrom)
  df$name <- str_c(df$chrom, df$summit, sep=":")
  return(df)})

# Assign 'high' or 'low' group label to each enhancer based on cluster
for(marker in MARKERS){
  enh_all[[marker]]$group <- ifelse(enh_all[[marker]]$cluster %in% clust_all$high[[marker]], "high", "low")
}

# Filter SSMs to keep only single nucleotide variants (SNVs)
SNVs <- SSMs %>% rename(., c("chromosome_start" = "start", "chromosome_end" = "end")) %>% 
  mutate(., chromosome = paste0("chr", .$chromosome)) %>%
  dplyr::filter(., end - start == 0)

# Filter SSMs to keep only single nucleotide variants (SNVs)
table_overlaps_markers <- list(
  "CtIP" = read_tsv(path_overlaps_ctip) %>% mutate(., group = ifelse(cluster %in% clust_all$high$CtIP, "high", "low")), 
  "GRHL" = read_tsv(path_overlaps_grhl) %>% mutate(., group = ifelse(cluster %in% clust_all$high$GRHL, "high", "low"))
) %>% suppressMessages()


##


# Find overlaps with random regions 

# Initialize lists to store random sequences and their overlaps with SNVs
input_variants <- SNVs
ran_seqs_df <- list()
ran_seqs_overlaps <- list()

# Loop over each marker (CtIP and GRHL)
for(marker in MARKERS){
  print(marker)
  marker_enh <- enh_all[[marker]]
  
  # Convert SNVs to GenomicRanges object
  SNVs_gr <- makeGRangesFromDataFrame(input_variants, keep.extra.columns = T)
  
  # Generate random sequences of the same length as the enhancer regions
  ran_seqs <- generate_random_seqs(SEED, dim(marker_enh)[1], win = WIN)
  ran_seqs$summit <- ran_seqs$start + (WIN/2)
  ran_seqs_df[[marker]] <- ran_seqs
  
  # Convert random sequences to GenomicRanges object
  ran_seqs_gr = makeGRangesFromDataFrame(ran_seqs, keep.extra.columns = T)
  
  # Find overlaps between random sequences and SNVs
  hits_obj_ran <- findOverlaps(query=ran_seqs_gr, subject=SNVs_gr)
  SNVs_sbj <- data.frame(SNVs_gr[subjectHits(hits_obj_ran)])
  colnames(SNVs_sbj)[1:5] <- paste(colnames(SNVs_sbj)[1:5], "sbj", sep = "_")
  ran_seqs_SNVs <- cbind(data.frame(ran_seqs_gr[queryHits(hits_obj_ran)], SNVs_sbj))
  
  # Store overlaps in the list
  ran_seqs_overlaps[[marker]] <- ran_seqs_SNVs
  
}


##


## Compute stats for multiple windows around summits

# Define different window sizes to analyze
WINS <- c(3000, 2000, 1000, 500, 100, 50)

# Initialize lists to store statistics for different windows
stat_enh_all_wins <- list()
stat_vars_all_wins <- list()

suppressMessages({

# Loop over each marker and window size
for(marker in MARKERS){
  
  for(win in WINS){
    
    # Extract overlaps for the current marker
    overlaps <- table_overlaps_markers[[marker]]
    input_enh <- enh_all[[marker]]
    
    # Filter overlaps to keep only those within the specified window
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
    
    
    ## Compute stats for random sequences
    
    # Extract overlaps and input random sequences data for the current marker
    overlaps <- ran_seqs_overlaps[[marker]]
    input_enh <- ran_seqs_df[[marker]]
    
    # Extract overlaps respecting window size 
    overlaps$dist <- abs(overlaps$start_sbj - overlaps$summit)  # Calculate distance between random sequence start and enhancer summit
    overlaps <- overlaps %>% dplyr::filter(., dist <= win)  # Filter for overlaps within the window size
    
    # Calculate the number of donors in which each enhancer is hit
    stat_enh_ran <- overlaps %>% dplyr::group_by(name) %>%
      dplyr::summarise(., n_donors_hit = length(unique(icgc_donor_id)))

    # Add to df enhancers with no mutations in any patient to the stats
    stat_enh_ran <- rbind(stat_enh_ran, data.frame("name" = input_enh$name[!input_enh$name %in% stat_enh_ran$name], 
                                               "n_donors_hit" = rep(0, length(input_enh$name[!input_enh$name %in% stat_enh_ran$name]))))
    
    # Calculate the number of variants within each enhancer region
    stat_vars_ran <- overlaps %>% group_by(name) %>% 
      dplyr::summarise(., n_vars = n()) # Count variants within each enhancer
    stat_vars_ran <- rbind(stat_vars_ran, data.frame("name" = input_enh$name[!input_enh$name %in% stat_vars_ran$name], 
                                             "n_vars" = rep(0, length(input_enh$name[!input_enh$name %in% stat_vars_ran$name]))))
    
    # Merge stats and save results
    stat_enh_all <- cbind.fill(stat_enh, stat_enh_ran) # Combine existing enhancer stats with random sequence stats
    colnames(stat_enh_all) <- c("name", "group", "n_donors_hit", "name.ran", "n_donors_hit_ran") # Rename columns
    stat_enh_all$n_donors_hit <- as.numeric(stat_enh_all$n_donors_hit)  # Ensure numeric type
    stat_enh_all$n_donors_hit_ran <- as.numeric(stat_enh_all$n_donors_hit_ran)  # Ensure numeric type
    
    stat_vars_all <- cbind.fill(stat_vars, stat_vars_ran) # Combine existing variant stats with random sequence stats
    colnames(stat_vars_all) <- c("name", "group", "n_vars", "name.ran", "n_vars_ran") # Rename columns
    stat_vars_all$n_vars <- as.numeric(stat_vars_all$n_vars) # Ensure numeric type
    stat_vars_all$n_vars_ran <- as.numeric(stat_vars_all$n_vars_ran) # Ensure numeric type
    
    # Save stats for the current window size and marker
    win_c <- as.character(win)
    stat_enh_all_wins[[marker]][[win_c]] <- stat_enh_all
    stat_vars_all_wins[[marker]][[win_c]] <- stat_vars_all
    
  }
}

})


##
  

suppressMessages({

## Plot stats
  
# Iterate over each marker
for(marker in MARKERS){
  # Iterate over each window size
  for(win in WINS){
    win_c <- as.character(win)  # Convert window size to character for naming
    stat_enh_to_plot <- stat_enh_all_wins[[marker]][[win_c]]  # Get enhancer stats for plotting
    stat_vars_to_plot <- stat_vars_all_wins[[marker]][[win_c]]  # Get variant stats for plotting
    
    # Compute the mean number of donors for each group
    avg <- stat_enh_to_plot %>% group_by(., group) %>% dplyr::summarise(., avg = mean(n_donors_hit))
    avg <- str_flatten(paste0(avg$group, ": ", round(avg$avg,3)), collapse = " ")
    
    # Compute the mean number of donors hit, comparing enhancers and random sequences
    mean_all <- stat_enh_to_plot %>% dplyr::select(., n_donors_hit, n_donors_hit_ran, group) %>%
      pivot_longer(-group) %>% group_by(., group, name) %>%
      dplyr::summarise(., mean = mean(value))
    
    # Plot number of donors hit per enhancer
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
    
    # Compute the mean number of variants per enhancer
    avg <- stat_vars_to_plot %>% group_by(., group) %>% dplyr::summarise(., avg = mean(n_vars))
    avg <- str_flatten(paste0(avg$group, ": ", round(avg$avg,3)), collapse = " ")
    
    # Compute the mean number of variants, comparing enhancers and random sequences
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
                         comparisons = list( c("high", "low")))+
      stat_compare_means(method = "wilcox", label = "p.signif")
    print(p2)
  }
}
  
})


##


