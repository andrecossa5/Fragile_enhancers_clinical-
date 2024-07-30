
library(tidyverse)
library(GenomicRanges)
#install.packages("/Users/ieo6983/Downloads/rowr_1.0.1.tar.gz", repos = NULL, type="source")
library(rowr)
library(ggpubr)

SEED <- 4321
set.seed(SEED)

source("/Users/ieo6983/Desktop/fragile_enhancer_clinical/utils/functions_genomics.R")

path_enhancers_ctip <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/downstream/CtIP_enh.hq_signal.clustered.tsv")
path_enhancers_grhl <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/downstream/GRHL_enh.hq_signal.clustered.tsv")

path_SSMs <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/genomics/pre_processed_ICGC/simple_somatic_mutation.open.matching_calls.with_AFs.tsv")
path_chrom_sizes <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/hg19.chrom.txt")

path_output_temp <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/results/ICGC/NEW/temp/")


MARKERS <- c("CtIP", "GRHL")
WIN <- 100
Mb <- 1000000

high_low <- T
save_temp <- F


##

# Clusters are the same for CtIP and GRHL
clust_high_ctip <- c("cluster_1", "cluster_2", "cluster_3", "cluster_4")
clust_low_ctip <- c("cluster_5")

clust_high <- list("CtIP" = clust_high_ctip, "GRHL" = clust_high_ctip)
clust_low <- list("CtIP" = clust_low_ctip, "GRHL" = clust_low_ctip)
clust_all <- list("high" = clust_high, "low" = clust_low)

##


# Define groups of clusters
cluster_groups <- setNames(vector(mode="list", length=length(MARKERS)), MARKERS)
cluster_groups$CtIP <- list(
  "high" = c("cluster_1", "cluster_2", "cluster_3", "cluster_4"), 
  "low" = c("cluster_5")
)
cluster_groups$GRHL <- cluster_groups$CtIP

# Read enhancers 
enh_ctip <- read_tsv(path_enhancers_ctip) 
enh_grhl <- read_tsv(path_enhancers_grhl) 
enh_all <- list("CtIP" = enh_ctip, "GRHL" = enh_grhl)

# add summit & name
enh_all <- lapply(enh_all, function(df){
  df$chrom <- paste0("chr", df$chrom)
  df$name <- str_c(df$chrom, df$summit, sep=":")
  return(df)})
# Add high-low cluster info
for(marker in MARKERS){
  enh_all[[marker]]$group <- ifelse(enh_all[[marker]]$cluster %in% clust_all$high[[marker]], "high", "low")
}


# Read variants
SSMs <- read_tsv(path_SSMs)
colnames(SSMs)[c(7,8)] <- c("start", "end")  
SSMs$chromosome <- paste("chr", SSMs$chromosome, sep = "")
SSMs <- SSMs[SSMs$end - SSMs$start == 0, ]

# Add mutatation ID equal to hartwig mutation ID
SSMs$ID <- str_sub(SSMs$chromosome, start = 4) %>% 
  str_c(., SSMs$end, sep = ":") %>% 
  str_c(., SSMs$reference_genome_allele, sep = "_") %>%
  str_c(., SSMs$mutated_to_allele, sep = "/")

SSMs_gr <- makeGRangesFromDataFrame(SSMs, keep.extra.columns = T)


# Read hg19 chrom sizes
standard_chrom <- paste0("chr", c(seq(1:22), "X", "Y"))
hg19_chrom_size <- read_tsv(path_chrom_sizes, col_names = c("chrom", "end")) %>% 
  dplyr::filter(., chrom %in% standard_chrom) %>%
  mutate(., start = 1) %>% relocate(., start, .before = end) 
hg19_chrom_size_gr <- makeGRangesFromDataFrame(hg19_chrom_size)


##

sample_var <- "icgc_sample_id"
fc_all_samples_markers <- list() # FD distribution for all CtIP/GRHL enhancers
fc_ran_all_samples_markers <- list() # FC distribution for random regions 
fc_x_group_all_samples_markers <- list() # FC distribution for high/low CtIP/GRHL2 enhancers 
back_mut_rate_all_samples <- c()

for(marker in MARKERS){
  
  # Select enhancers 
  marker_enh <- enh_all[[marker]]
  marker_enh$start <- marker_enh$summit - WIN
  marker_enh$end <- marker_enh$summit + WIN
  marker_enh$length <- marker_enh$end - marker_enh$start
  enh_gr <- makeGRangesFromDataFrame(marker_enh, keep.extra.columns = T)
  
  # Store FoldChange(mutation rate over background)
  fc_all_samples <- c()
  fc_ran_all_samples <- c()
  fc_x_group_all_samples <- list("high" = c(), "low" = c())
  
  for(i in 1:length(unique(SSMs[[sample_var]]))){
    sample <- unique(SSMs[[sample_var]])[i]
    print(paste0("Analyzing sample: ", sample))
  
    # Generate random genomic sequences
    ran_seqs <- generate_random_seqs(SEED+i, dim(marker_enh)[1], win = WIN)
    ran_seqs_gr = makeGRangesFromDataFrame(ran_seqs, keep.extra.columns = T)
    
    # Get SNVs from 1 sample only 
    sample_SSMs_gr <- SSMs_gr[mcols(SSMs_gr)[[sample_var]] == sample]
    
    #hits_background
    hits_background <- findOverlaps(query=hg19_chrom_size_gr, subject=sample_SSMs_gr)
    sample_back_SSMs <- cbind(data.frame(hg19_chrom_size_gr[queryHits(hits_background)]), data.frame(sample_SSMs_gr[subjectHits(hits_background)]))
    tot_back_variants <- length(unique(sample_back_SSMs$ID))
    tot_back_length_Mb <- sum(hg19_chrom_size$end) / Mb
    # Compute background mutation frequency x Mb: number of SNVs x Mb (whole genome)     
    back_mut_rate <- tot_back_variants / tot_back_length_Mb
    back_mut_rate_all_samples <- c(back_mut_rate_all_samples, back_mut_rate)
    
    ##hits_enh
    hits_enh <- findOverlaps(query=enh_gr, subject=sample_SSMs_gr)
    print(paste0("Number of hits found for enhancers: ", length(hits_enh)))
    sample_enh_SSMs <- cbind(data.frame(enh_gr[queryHits(hits_enh)]), data.frame(sample_SSMs_gr[subjectHits(hits_enh)]))
    tot_enh_variants <- length(unique(sample_enh_SSMs$ID))
    tot_enh_length_Mb <- sum(marker_enh$end - marker_enh$start) / Mb
    # Compute enhancer mutation frequency x Mb: number of SNVs x Mb (enhancers)     
    enh_mut_rate <- tot_enh_variants / tot_enh_length_Mb
    
    # X group (high - low) check mutation rate x Mb over background 
    if(high_low == T){
      if(!dim(sample_enh_SSMs)[1] == 0){
        tot_enh_variants <- sample_enh_SSMs %>% as_tibble(.name_repair = "unique") %>%
          group_by(., group) %>%
          dplyr::summarise(., tot_enh_var_x_clust = length(unique(ID)))
        # total length of high-union and low-union    
        length_Mb <- marker_enh %>% group_by(., group) %>% 
          dplyr::summarise(., tot_len = sum(length), .groups = "keep") %>%
          mutate(., tot_len_x_clust_Mb = tot_len / Mb)
        
        # Compute mut_rate for each cluster
        mut_rate <- left_join(tot_enh_variants, length_Mb, by = "group") %>% mutate(., enh_mut_rate = tot_enh_var_x_clust / tot_len_x_clust_Mb) %>% 
          arrange(., group)
        
        # Fold change over background
        mut_rate_high <- mut_rate[mut_rate$group == "high", ]$enh_mut_rate / back_mut_rate
        mut_rate_high <- ifelse(length(mut_rate_high) == 0, 0, mut_rate_high)
        mut_rate_low <- mut_rate[mut_rate$group == "low", ]$enh_mut_rate / back_mut_rate
        mut_rate_low <- ifelse(length(mut_rate_low) == 0, 0, mut_rate_low)
        
        fc_x_group_all_samples[["high"]] <- c(fc_x_group_all_samples[["high"]], mut_rate_high)
        fc_x_group_all_samples[["low"]] <- c(fc_x_group_all_samples[["low"]], mut_rate_low)
        
      } else {
        tot_enh_variants <- c(0,0)
        mut_rate <- c(0,0)
        
        # Fold change over background
        fc_x_group_all_samples[["high"]] <- c(fc_x_group_all_samples[["high"]], 0)
        fc_x_group_all_samples[["low"]] <- c(fc_x_group_all_samples[["low"]], 0)
      }
    }
    
    # Fold change ( enh_mut_rate / back_mut_rate)
    fc <- enh_mut_rate / back_mut_rate
    fc_all_samples <- c(fc_all_samples, fc)
    
    ##hits random 
    hits_ran <- findOverlaps(query=ran_seqs_gr, subject=sample_SSMs_gr)
    print(paste0("Number of hits found for random seqs: ", length(hits_ran)))
    ran_SSMs <- cbind(data.frame(ran_seqs_gr[queryHits(hits_ran)]), data.frame(SSMs_gr[subjectHits(hits_ran)]))
    tot_ran_variants <- length(unique(ran_SSMs$ID))
    tot_ran_length_Mb <- sum(ran_seqs$end - ran_seqs$start) / Mb
    # Compute enhancer mutation frequency x Mb: number of SNVs x Mb (enhancers)     
    ran_mut_rate <- tot_ran_variants / tot_ran_length_Mb
    
    # Fold change ( ran_mut_rate / back_mut_rate)
    fc_ran <- ran_mut_rate / back_mut_rate
    fc_ran_all_samples <- c(fc_ran_all_samples, fc_ran)
    
  }
  fc_all_samples_markers[[marker]] <- fc_all_samples
  fc_ran_all_samples_markers[[marker]] <- fc_ran_all_samples
  fc_x_group_all_samples_markers[[marker]] <- fc_x_group_all_samples
}


##


# Plot FC distribution for all samples
for(marker in MARKERS){
  
  # Enhancers mutation rate only
  df_to_plot <- data.frame("fc_enh" = fc_all_samples_markers[[marker]], 
                           "fc_ran" = fc_ran_all_samples_markers[[marker]])
  p <- df_to_plot %>% pivot_longer(everything(), names_to = "seq_type", values_to = "mut_rate") %>%
    ggplot(., aes(x = seq_type, y = mut_rate, group = seq_type))+
    geom_jitter(width = 0.3, alpha = 0.4)+
    geom_boxplot(alpha=0.2, width = 0.4, outliers = F)+
    stat_compare_means(label = "p.signif", method = "wilcox.test", size = 3, comparisons = list(c("fc_enh", "fc_ran")))+
    theme_light()+
    geom_hline(yintercept = 1, col = "grey", linetype = 2)+
    labs(x = "", y = "Fold-change over background", 
         title = paste0("Mutation frequency in ", marker, " enhancers"))
  print(p)
  # Same but capped 
  p <- df_to_plot %>% pivot_longer(everything(), names_to = "seq_type", values_to = "mut_rate") %>%
    ggplot(., aes(x = seq_type, y = mut_rate, group = seq_type))+
    geom_jitter(width = 0.3, alpha = 0.4)+
    geom_boxplot(alpha=0.2, width = 0.4, outliers = F)+
    stat_compare_means(label = "p.signif", method = "wilcox.test", size = 3, comparisons = list(c("fc_enh", "fc_ran")))+
    theme_light()+
    geom_hline(yintercept = 1, col = "grey", linetype = 2)+
    labs(x = "", y = "Fold-change over background", 
         title = paste0("Mutation frequency in ", marker, " enhancers"))+
    ylim(c(0,10))
  print(p)
  
  # Same jitter_plot, but in log2 scale
  df_to_plot_log <- data.frame("fc_enh" = log2(fc_all_samples_markers[[marker]]+0.1), 
                           "fc_ran" = log2(fc_ran_all_samples_markers[[marker]]+0.1))
  p <- df_to_plot_log %>% pivot_longer(everything(), names_to = "seq_type", values_to = "mut_rate") %>%
    ggplot(., aes(x = seq_type, y = mut_rate, group = seq_type))+
    geom_jitter(width = 0.3, alpha = 0.4)+
    geom_boxplot(alpha=0.2, width = 0.4, outliers = F)+
    stat_compare_means(label = "p.signif", method = "wilcox.test", size = 3, comparisons = list(c("fc_enh", "fc_ran")))+
    theme_light()+
    geom_hline(yintercept = 0, col = "grey", linetype = 2)+
    labs(x = "", y = "log2(Fold-change over background)", 
         title = paste0("Mutation frequency in ", marker, " enhancers"))
  print(p)
  
  
  ## 
  
  
  # Per group (high-low)
  rownames(fc_all_samples_markers[[marker]]) <- NULL
  df_to_plot <- cbind.fill(fc_all_samples_markers[[marker]], 
             fc_ran_all_samples_markers[[marker]], 
             fc_x_group_all_samples_markers[[marker]][["high"]], 
             fc_x_group_all_samples_markers[[marker]][["low"]],
             fill = NA)
  colnames(df_to_plot) <- c("fc_enh", "fc_ran", "fc_high", "fc_low")
  
  p <- df_to_plot %>% pivot_longer(everything(), names_to = "seq_type", values_to = "mut_rate") %>%
    ggplot(., aes(x = seq_type, y = mut_rate, group = seq_type))+
    geom_jitter(width = 0.3, alpha = 0.4)+
    geom_boxplot(alpha=0.2, width = 0.4, outliers = F)+
    stat_compare_means(label = "p.signif", method = "wilcox.test", size = 3, 
                       comparisons = list(c("fc_enh", "fc_ran"), c("fc_high","fc_ran"), c("fc_low", "fc_ran"), c("fc_high", "fc_low")))+
    theme_light()+
    geom_hline(yintercept = 1, col = "grey", linetype = 2)+
    labs(x = "", y = "Fold-change over background", 
         title = paste0("Mutation frequency in ", marker, " enhancers"))+
    stat_summary(
      fun = mean, 
      geom = "text", 
      aes(label = sprintf("%.2f", ..y..)), 
      vjust = -20,
      color = "blue"
    )
  print(p)
  
  p <- df_to_plot %>% pivot_longer(everything(), names_to = "seq_type", values_to = "mut_rate") %>%
    ggplot(., aes(x = seq_type, y = mut_rate, group = seq_type))+
    geom_jitter(width = 0.3, alpha = 0.4)+
    geom_boxplot(alpha=0.2, width = 0.4, outliers = F)+
    theme_light()+
    geom_hline(yintercept = 1, col = "grey", linetype = 2)+
    labs(x = "", y = "Fold-change over background", 
         title = paste0("Mutation frequency in ", marker, " enhancers"))+
    ylim(c(0,15))
  print(p)
  
  # Same jitter_plot, but in log2 scale
  df_to_plot_log <- log2(df_to_plot+0.1)
  p <- df_to_plot_log %>% pivot_longer(everything(), names_to = "seq_type", values_to = "mut_rate") %>%
    ggplot(., aes(x = seq_type, y = mut_rate, group = seq_type))+
    geom_jitter(width = 0.3, alpha = 0.4)+
    geom_boxplot(alpha=0.2, width = 0.4, outliers = F)+
    stat_compare_means(label = "p.signif", method = "wilcox.test", size = 3, 
                       comparisons = list(c("fc_enh", "fc_ran"), c("fc_high","fc_ran"), c("fc_low", "fc_ran"), c("fc_high", "fc_low")))+
    theme_light()+
    geom_hline(yintercept = 1, col = "grey", linetype = 2)+
    labs(x = "", y = "log2(Fold-change over background)", 
         title = paste0("Mutation frequency in ", marker, " enhancers"))+
    stat_summary(
      fun = mean, 
      geom = "text", 
      aes(label = sprintf("%.2f", ..y..)), 
      vjust = -18.5,
      color = "blue"
    )
  print(p)
  
  # Same jitter_plot, in log2 scale, capped 
  df_to_plot_log <- log2(df_to_plot+0.1)
  p <- df_to_plot_log %>% pivot_longer(everything(), names_to = "seq_type", values_to = "mut_rate") %>%
    ggplot(., aes(x = seq_type, y = mut_rate, group = seq_type))+
    geom_jitter(width = 0.3, alpha = 0.4)+
    geom_boxplot(alpha=0.2, width = 0.4, outliers = F)+
    theme_light()+
    geom_hline(yintercept = 1, col = "grey", linetype = 2)+
    labs(x = "", y = "log2(Fold-change over background)", 
         title = paste0("Mutation frequency in ", marker, " enhancers"))+
    ylim(c(-2,13))+
    stat_summary(
      fun = mean, 
      geom = "text", 
      aes(label = sprintf("%.2f", ..y..)), 
      vjust = -25,
      color = "blue"
    )+
    stat_compare_means(label = "p.signif", method = "wilcox.test", size = 3, 
                       comparisons = list(c("fc_enh", "fc_ran"), c("fc_high","fc_ran"), c("fc_low", "fc_ran"), c("fc_high", "fc_low")))
  print(p)
}


##


if(save_temp == T){
  # Save ICGC objects temporarily  
  saveRDS(fc_all_samples_markers, file = fs::path(path_output_temp, "fc_dist.all_enhancers.icgc.rds"))
  saveRDS(fc_ran_all_samples_markers, file = fs::path(path_output_temp, "fc_dist.all_random.icgc.rds"))
  saveRDS(fc_x_group_all_samples_markers, file = fs::path(path_output_temp, "fc_dist.grouped_enhancers.icgc.rds"))
}


##


# TODO: 

# - implement for hartwig as well
# - implement for hartwig vs. ICGC

# - Exclude hypermutated regions under physiological conditions?

