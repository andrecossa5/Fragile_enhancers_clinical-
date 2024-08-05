
### REMOVE_SELF_OVERLAPS
# Given a GRages object containing a set of genomic regions, remove overlapping regions within the GRanges 

remove_self_overlaps <- function(input_granges){
  
  self_ovrlps <- data.frame(findOverlaps(input_granges)) # overlap GRanges_obj with itself
  self_ovrlps$hit <- self_ovrlps$queryHits == self_ovrlps$subjectHits 
  
  filter_overlaps <- self_ovrlps %>% dplyr::group_by(queryHits) %>%
    dplyr::arrange(hit, .by_group = T) %>% # arrange by T/F
    filter(row_number() == 1) %>% # and filter out regions overlapping with themselves only
    ungroup() %>% dplyr::select(hit)
  output_granges <- input_granges[filter_overlaps$hit] 
  
  return(output_granges)
}





### GENERATE RANDOM SEQS ###
# Generate a set of random genomic sequences 

generate_random_seqs <- function(seed, n_seqs, win, chrom_lengths_info){
  set.seed(seed)
  library(tidyverse)
  path_chrom_sizes <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/hg19.chrom.txt")
  
  if(file.exists(path_chrom_sizes) == F){
    print("File with chromosome sizes does not exist, using hardcoded sizes:"); cat("\n")
    
    standard_chrom <- paste0("chr", c(seq(1:22), "X", "Y"))
    seq_lengths <- c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 155270560, 146364022, 141213431, 135534747, 
                     135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 63025520, 59373566, 59128983, 51304566, 48129895)
    chrs_info <- data.frame("chrom" = standard_chrom, "seq_lengths" = seq_lengths)
  } else {
    # Define chromosomes info
    standard_chrom <- paste0("chr", c(seq(1:22), "X", "Y"))
    chrs_info <- read_tsv(path_chrom_sizes, col_names = c("chrom", "seq_lengths")) %>% 
      dplyr::filter(., chrom %in% standard_chrom) %>% suppressMessages()
  }
  
  # Sample chromosomes and start positions 
  ran_seqs <- data.frame(chrom = sample(chrs_info$chrom, n_seqs, replace=T))
  start_pos <- c()
  for(i in 1:length(ran_seqs$chrom)){
    chrom <- ran_seqs$chrom[i]
    chrom_length <- chrs_info[chrs_info$chrom == chrom, ]$seq_lengths
    start_pos[i] <- round(runif(1, 0, chrom_length),0)
  }
  ran_seqs$start <- start_pos
  
  # Add end position and sequence name
  ran_seqs$end <- ran_seqs$start + (win*2)
  ran_seqs$name <- paste(ran_seqs$chrom, ran_seqs$start, ran_seqs$end, sep = "_")
  
  return(ran_seqs)  
}




### STAT ENHANCERS ###
# Computes number of donors hit by each enhancer 
# Inputs must be 2 lists: 
# - a list of data.frames with overlaps 
# - a list of the dataframes containing the data.frames containing the enhancers used for the overlaps
# classes = vector containing classes of enhancers/sequences overlapped

compute_stat_enhancers <- function(input_overlaps, source_enhancers_df, classes){
  df_all_stats <- data.frame(matrix(ncol = length(input_overlaps), # Bottom line will give error sapply(list(source_enhancers_df)) ??
                                    nrow = max(sapply(source_enhancers_df, function(df) dim(df)[1]))))
  
  colnames(df_all_stats) <- classes
  
  for(i in 1:length(input_overlaps)){
    df_overlaps <- input_overlaps[[i]]
    df_source <- source_enhancers_df[[i]]
    
    # Compute number of donors in which each enhancer is mutated
    stat_enh <- df_overlaps %>% dplyr::group_by(name) %>%
      dplyr::summarise(., n_donors_hit = length(unique(icgc_donor_id)))
    
    # Add to df enhancers with no mutation in any patient
    stat_enh <- rbind(stat_enh, 
                      data.frame("name" = df_source$name[!df_source$name %in% stat_enh$name], 
                                 "n_donors_hit" = rep(0, length(df_source$name[!df_source$name %in% stat_enh$name]))))
    
    # add to df_all_stats
    df_all_stats[, i] <- c(stat_enh$n_donors_hit, 
                           rep(NA, length(df_all_stats[, i]) - length(stat_enh$n_donors_hit)))
  }
  
  return(df_all_stats)
}





### STAT DONORS
# Compute number of enhancers hit/mutated in each donor
# input_overlaps = df of overlaps
# source_donors_df = df with all donors used up to now

compute_stat_donors <- function(input_overlaps, source_donors_df, classes){
  df_all_stats <- data.frame(matrix(ncol = length(input_overlaps), 
                                    nrow = length(unique(source_donors_df$icgc_donor_id))))
  colnames(df_all_stats) <- classes
  
  for(i in 1:length(input_overlaps)){
    df_overlaps <- input_overlaps[[i]]
    
    # Compute number of enhancers mutated in each donor
    stat_donors <- df_overlaps %>% group_by(icgc_donor_id) %>%
      dplyr::summarise("n_enhancers_hit" = length(unique(name)))
    
    # Add to df donors in which no enhancer is mutated
    stat_donors <- rbind(stat_donors, 
                         data.frame("icgc_donor_id" = unique(source_donors_df$icgc_donor_id[!source_donors_df$icgc_donor_id %in% stat_donors$icgc_donor_id]), 
                                    "n_enhancers_hit" = rep(0, length(unique(source_donors_df$icgc_donor_id[!source_donors_df$icgc_donor_id %in% stat_donors$icgc_donor_id])))))
    # add to df_all_stats
    df_all_stats[, i] <- stat_donors$n_enhancers_hit
  }
  
  return(df_all_stats)
}



### Compute P-values distribution for case-control comparison ###
# n_samples = numero di mutliple samplings
# start_seed = starting seed. For each i in n_samples draws, start_seed+i will be used to sample
# cluster_enh = a dataframe containing the group of enhancers, belonging to a certain cluster, representing the 'case' group. 
  # the dim()[1] of the dataframe is used to sample the right number of random sequences.
# input_variants_gr = GRanges object with sequences of the variants of interest, for which overlaps with random regions will be searched.
# n_mutations_case = number of mutations already found for the 'case' group of enahncers, to be compared with each random sample of random regions
# (should generalize this one to be able to use different variables, not only n_muts)

compute_pval_dist_ctrl_random <- function(n_samples = 1000, start_seed = 4321, 
                                          cluster_enh, input_variants_gr, 
                                          n_mutations_case){
  # i.e. cluster = marker_enh_clust
  p_val_dist <- c()
  
  for(i in 1:n_samples){
    ran_seqs <- generate_random_seqs(seed=start_seed+i, n_seqs=dim(cluster_enh)[1], win = WIN)
    ran_seqs_gr = makeGRangesFromDataFrame(ran_seqs, keep.extra.columns = T)
    
    # Random
    hits_obj_ran <- findOverlaps(query=ran_seqs_gr, subject=input_variants_gr)
    SSMs_sbj <- data.frame(input_variants_gr[subjectHits(hits_obj_ran)])
    colnames(SSMs_sbj)[1:5] <- paste(colnames(SSMs_sbj)[1:5], "sbj", sep = "_")
    ran_seqs_SSMs <- cbind(data.frame(ran_seqs_gr[queryHits(hits_obj_ran)], SSMs_sbj))
    
    # Count number of mutations overlapping with each random seq
    n_muts_ran <- ran_seqs_SSMs %>% group_by(name) %>% 
      dplyr::summarise(n_muts = n()) %>%
      ungroup() %>%
      dplyr::select(name, n_muts) 
    # Add enhancers with 0 mutations
    ran_0_muts <- unique(ran_seqs$name[!ran_seqs$name %in% n_muts_ran$name])
    n_muts_ran <- rbind(n_muts_ran, data.frame("name" = ran_0_muts, "n_muts" = 0))
    n_muts_ran <- arrange(n_muts_ran, desc(n_muts_ran$n_muts))
    
    p_val <- round(wilcox.test(n_mutations_case$n_muts, n_muts_ran$n_muts)$p.value,3)
    p_val_dist <- c(p_val_dist, p_val)
    
  }
  
  return(p_val_dist)
}




### P-VALUES CASES ###
compute_pval_dist_cases <- function(n_samples = 1000, start_seed = 4321, 
                                   sample_size = 1900, # number of enhancers to resample from each group
                                   df_vars){ # dataframe with variable to compare for each group, i.e. 'n_muts
  p_val_dist <- list()
  for(marker in MARKERS){
    p_val_dist[[marker]] <- c()
    
    for(i in 1:n_samples){
      set.seed(start_seed+i)
      
      # Extract variable of interest for each group
      var_high <- na.omit(df_vars[[marker]] %>% dplyr::select("high"))
      var_low <- na.omit(df_vars[[marker]] %>% dplyr::select("low"))
      
      # Sample X observations with replacement 
      sample_high <- sample(x = 1:dim(var_high)[1], size = sample_size, replace = T)
      sample_high <- var_high[sample_high, ]
      sample_low <- sample(x = 1:dim(var_low)[1], size = sample_size, replace = T)
      sample_low <- var_low[sample_low, ]
      
      # Compute significance of difference in 'var' (i.e. number of mutations) among the 2 groups      
      p_val <- round(wilcox.test(sample_high, sample_low)$p.value,3)
      p_val_dist[[marker]] <- c(p_val_dist[[marker]], p_val)
    } 
  }
  
  # Return df with p-values for high-low comparison, for each marker 
  return(p_val_dist)
} 



#' Plot distribution of SNVs over enhancers
#' 
#' Given ...
#' 
#' @param 
#' @returns description
#' TODO: implement 


