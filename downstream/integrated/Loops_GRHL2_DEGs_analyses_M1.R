
library(tidyverse)
library(GenomicRanges)
library(viridis)
library(gridExtra)


### Hi-ChIP Loops ###
kb <- 2
location <- "local" # 'local' or 'hpc'

# Source custom functions for processing Hi-ChIP loops
source("/Users/ieo6983/Desktop/fragile_enhancer_clinical/utils/loops_functions.R")

# Define file paths for data based on location
path_hichip <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/HiChip/filtered_loops/")
path_enhancers <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/downstream/")
path_enhancers_ctip <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/downstream/CtIP_enh.hq_signal.clustered.tsv")
path_enhancers_grhl <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/downstream/GRHL_enh.hq_signal.clustered.tsv")
path_tss <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/others/TSSs_elisa/TSSs_from_USCS_hg19_EMSEMBL.tsv")
path_degs <- fs::path("/Users/ieo6983/Desktop/expression/DEGs_length_scaled/Df_DEGs.df_LFC_sig.padj_0.05.log2FC_1.Up_and_Down.tsv")

# MSigDB Gene Set - Can be downloaded at: https://www.gsea-msigdb.org/gsea/msigdb/human/genesets.jsp
path_hallmark <- fs::path("/Users/ieo6983/Desktop/expression/GSEA/DB/hallmark_gene_sets.h.all.v2023.2.Hs.symbols.gmt")

# Define output paths
path_anno_loops <- fs::path(paste0("/Users/ieo6983/Desktop/fragile_enhancer_clinical/results/integrated/NEW/", kb, "kb/data/"))
path_results <- fs::path(paste0("/Users/ieo6983/Desktop/fragile_enhancer_clinical/results/integrated/NEW/", kb, "kb"))

# Adjust paths if running on HPC
if(location == "hpc"){
  source("/hpcnfs/scratch/PGP/Ciacci_et_al/fragile_enhancer_clinical/utils/loops_functions.R")
  path_hichip <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/HiChip/filtered_loops/")
  path_enhancers <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/downstream/")
  path_enhancers_ctip <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/downstream/CtIP_enh.hq_signal.clustered.tsv")
  path_enhancers_grhl <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/downstream/GRHL_enh.hq_signal.clustered.tsv")
  
  path_tss <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/others/TSSs_elisa/TSSs_from_USCS_hg19_EMSEMBL.tsv")
  path_degs <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/results/expression/RNA_seq/DEGs/Df_DEGs.df_LFC_sig.padj_0.05.log2FC_1.Up_and_Down.tsv")
  path_hallmark <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/results/integrated/NEW/hallmark_gene_sets.h.all.v2023.2.Hs.symbols.gmt")
  
  path_anno_loops <- fs::path(paste0("/hpcnfs/scratch/PGP/Ciacci_et_al/results/integrated/NEW/", kb, "kb/data/"))
  path_results <- fs::path(paste0("/hpcnfs/scratch/PGP/Ciacci_et_al/results/integrated/NEW/", kb, "kb"))
}


##


# Read loop files for enhancers associated with DEGs for SCR and KD conditions
scr_enh_degs_only <- read_tsv(fs::path(path_anno_loops, sprintf("%skb_SCR.anno_loops.GRHL2_enh_DEGs_prom_ONLY.tsv", kb)))
kd_enh_degs_only <- read_tsv(fs::path(path_anno_loops, sprintf("%skb_KD.anno_loops.GRHL2_enh_DEGs_prom_ONLY.tsv", kb)))

# Read condition-specific loops for SCR and KD
scr_specific <- read_tsv(fs::path(path_anno_loops, sprintf("%skb_SCR_specific_loops.tsv", kb)))
kd_specific <- read_tsv(fs::path(path_anno_loops, sprintf("%skb_KD_specific_loops.tsv", kb)))

# Read GRHL2-bound enhancers
enh_grhl2 <- read_tsv(path_enhancers_grhl)

# Update enhancer names to include chromosome prefix and summit information
enh_grhl2$chrom <- paste0("chr", enh_grhl2$chrom)
enh_grhl2$name <- str_c(enh_grhl2$chrom, enh_grhl2$summit, sep=":")

# Categorize enhancers based on cluster assignment
clusters_high <- c("cluster_1", "cluster_2", "cluster_3", "cluster_4")
clusters_low <- c("cluster_5")
enh_grhl2$group <- "high"
enh_grhl2[enh_grhl2$cluster %in% clusters_low, ]$group <- "low"


##


# Unbundle ambiguous loops
scr_enh_degs_only_filt <- unbundle_ambiguous_bins(scr_enh_degs_only)
kd_enh_degs_only_filt <- unbundle_ambiguous_bins(kd_enh_degs_only)


##


# Create a list of filtered loop data for SCR and KD conditions
enh_degs_only <- list("scr" = scr_enh_degs_only_filt, "kd" = kd_enh_degs_only_filt)

# Count associations between DEGs and loops for SCR condition
t1 <- table(c(scr_enh_degs_only_filt$DE1, scr_enh_degs_only_filt$DE2)[!is.na(c(scr_enh_degs_only_filt$DE1, scr_enh_degs_only_filt$DE2))])

# Count associations between DEGs and loops for KD condition
t2 <- table(c(kd_enh_degs_only_filt$DE1, kd_enh_degs_only_filt$DE2)[!is.na(c(kd_enh_degs_only_filt$DE1, kd_enh_degs_only_filt$DE2))])

# Compile counts into a data frame for visualization
loop_degs <- data.frame(group = c("scr", "kd"),
                        down = c(t1["Down"], t2["Down"]),
                        up = c(t1["Up"], t2["Up"]))
loop_degs$group <- factor(loop_degs$group, levels = c("scr", "kd"))

# Plot the number of loops associated with upregulated and downregulated DEGs
loop_degs %>% pivot_longer(-group) %>%
  ggplot(., aes(x=group, y = value, fill = name))+
  geom_bar(stat = "identity", position = "dodge")+
  theme_light()+
  labs(title = "Loops association to DEGs", 
       x = "", y = "number of loops")+
  scale_fill_manual(values = c("purple3", "orange3"))+
  theme(axis.text.x = element_text(size = 12, color = "black"))+
  scale_y_continuous(n.breaks = 10)


# Read DEG data
DEGs <- read_tsv(path_degs)

# Initialize list to store DEGs involved in loops
degs_in_loops <- list()

# Loop through each condition (SCR and KD) to find DEGs involved in loops
for(cond in c("scr", "kd")){
  df_loops <- enh_degs_only[[cond]]
  
  # Extract unique DEGs from the loop data
  degs_uniq <- data.frame("gene_name" = unique(c(na.omit(df_loops$gene_name1), na.omit(df_loops$gene_name2))))
  degs_uniq <- left_join(degs_uniq, DEGs[, c("gene_name", "DE", "log2FoldChange", "padj")], by = "gene_name")
  
  degs_in_loops[[cond]] <- degs_uniq
  # Optionally save DEG data
  #degs_uniq %>% write_tsv(., fs::path(path_results, sprintf("/data/%skb_%s.DEGs_in_enh_deg_loops",kb, toupper(cond)), ext = "tsv"))
}

# Count the number of unique DEGs involved in loops for each condition
for(cond in c("scr", "kd")){
  dfs_cond <- enh_degs_only[[cond]]
  
  # Count unique DEGs in loops
  degs_n_uniq <- length(unique(c(na.omit(dfs_cond$gene_name1), na.omit(dfs_cond$gene_name2))))
  print(sprintf("%s: %s", cond, degs_n_uniq))
  
  # Count upregulated and downregulated DEGs in loops
  up1 <- dfs_cond %>% filter(DE1 == "Up")
  up2 <- dfs_cond %>% filter(DE2 == "Up")
  degs_up <- length(unique(na.omit(c(up1$gene_name1, up2$gene_name2))))
  down1 <- dfs_cond %>% filter(DE1 == "Down")
  down2 <- dfs_cond %>% filter(DE2 == "Down")
  degs_down <- length(unique(na.omit(c(down1$gene_name1, down2$gene_name2))))
  print(sprintf("%s - Down: %s; Up: %s", cond, degs_down, degs_up))
}


# Analyze cluster composition of enhancers involved in loops
for(cond in c("scr", "kd")){
  df_loops <- enh_degs_only[[cond]]
  
  # Combine clusters from loop data and exclude NA values
  clusters_in_loops <- c(df_loops$cluster1, df_loops$cluster2)[!is.na(c(df_loops$cluster1, df_loops$cluster2))]
  
  # Print cluster composition
  print(sprintf("Cluster composition of ENH-DEG loops - %s", toupper(cond)))
  print(table(clusters_in_loops))
  cat("\n")
  
  # Calculate high and low cluster proportions
  h <- as.data.frame(table(clusters_in_loops)) %>% 
    dplyr::filter(clusters_in_loops %in% clusters_high) %>% 
    dplyr::select(Freq) %>% colSums()
  l <- as.data.frame(table(clusters_in_loops)) %>% 
    dplyr::filter(clusters_in_loops %in% clusters_low) %>% 
    dplyr::select(Freq) %>% colSums()
  
  hr <- h / sum(enh_grhl2$group == "high")
  lr <- l / sum(enh_grhl2$group == "low")
  
  # Print proportions and perform chi-square test
  print(sprintf("Number of enh. in clusters HIGH: %s (%s) vs. number of enh. in clusters LOW: %s (%s)", 
                h, round(hr,3), 
                l, round(lr,3)))
  cat("\n")
  
  # Define observed counts for chi-square test
  observed_counts <- matrix(c(h, l, sum(enh_grhl2$group == "high")-h, sum(enh_grhl2$group == "low")-l), nrow = 2, byrow = TRUE)
  colnames(observed_counts) <- c("High", "Low")
  rownames(observed_counts) <- c("In loops", "Not in loops")
  chi_squared_test <- chisq.test(observed_counts)
  
  # Print chi-square test p-value
  print(sprintf("Chi-square test p-value: %s", round(chi_squared_test$p.value,3)))
  cat("\n")
}


##


# SCR-specific enhancers connected to DEGs
# Identify and filter SCR-specific enhancers that are connected to DEGs
scr_enh_degs_only_spec <- inner_join(scr_enh_degs_only_filt, scr_specific, by = c('seqnames1', 'start1', 'end1', 'seqnames2', 'start2', 'end2'))

# Create a table of DEG categories (Up/Down) for SCR-specific enhancers
t1 <- table(c(scr_enh_degs_only_spec$DE1, scr_enh_degs_only_spec$DE2)[!is.na(c(scr_enh_degs_only_spec$DE1, scr_enh_degs_only_spec$DE2))])

# KD-specific enhancers connected to DEGs
# Identify and filter KD-specific enhancers that are connected to DEGs
kd_enh_degs_only_spec <- inner_join(kd_enh_degs_only_filt, kd_specific, by = c('seqnames1', 'start1', 'end1', 'seqnames2', 'start2', 'end2'))

# Create a table of DEG categories (Up/Down) for KD-specific enhancers
t2 <- table(c(kd_enh_degs_only_spec$DE1, kd_enh_degs_only_spec$DE2)[!is.na(c(kd_enh_degs_only_spec$DE1, kd_enh_degs_only_spec$DE2))])

# Combine the results into a data frame for plotting
cond_spec_degs <- data.frame(group = c("scr", "kd"),
                             down = c(t1["Down"], t2["Down"]), 
                             up = c(t1["Up"], t2["Up"]))
cond_spec_degs$group <- factor(cond_spec_degs$group, levels = c("scr", "kd"))

# Plot the number of condition-specific loops associated with upregulated and downregulated DEGs
cond_spec_degs %>% pivot_longer(-group) %>%
  ggplot(., aes(x=group, y = value, fill = name))+
  geom_bar(stat = "identity", position = "dodge")+
  theme_light()+
  labs(title = "Condition-specific loops & their association to DEGs", 
       x = "", y = "number of loops")+
  scale_fill_manual(values = c("purple3", "orange3"))+
  theme(axis.text.x = element_text(size = 12, color = "black"))

# Number of DEGs involved in condition-specific loops
cond_spec <-  list("scr" = scr_enh_degs_only_spec, "kd" = kd_enh_degs_only_spec)

# Iterate through conditions (SCR and KD) to calculate the number of unique DEGs
for(cond in c("scr", "kd")){
  dfs_cond <- cond_spec[[cond]]
  
  # Count unique DEGs involved in loops
  degs_n_uniq <- length(unique(c(na.omit(dfs_cond$gene_name1), na.omit(dfs_cond$gene_name2))))
  print(sprintf("%s: %s", cond, degs_n_uniq))
  
  # Count upregulated DEGs
  up1 <- dfs_cond %>% filter(DE1 == "Up")
  up2 <- dfs_cond %>% filter(DE2 == "Up")
  degs_up <- length(unique(na.omit(c(up1$gene_name1, up2$gene_name2))))
  
  # Count downregulated DEGs
  down1 <- dfs_cond %>% filter(DE1 == "Down")
  down2 <- dfs_cond %>% filter(DE2 == "Down")
  degs_down <- length(unique(na.omit(c(down1$gene_name1, down2$gene_name2))))
  
  # Print counts of upregulated and downregulated DEGs
  print(sprintf("%s - Down: %s; Up: %s", cond, degs_down, degs_up))
}

# DEGs involved in loops
DEGs <- read_tsv(path_degs)

# Initialize list to store DEGs in condition-specific loops
degs_in_loops_spec <- list()

# Iterate through conditions (SCR and KD) to find DEGs in condition-specific loops

for(cond in c("scr", "kd")){
  df_loops <- cond_spec[[cond]]
  
  # Create a data frame with unique DEGs and their associated attributes
  degs_uniq <- data.frame("gene_name" = unique(c(na.omit(df_loops$gene_name1), na.omit(df_loops$gene_name2))))
  degs_uniq <- left_join(degs_uniq, DEGs[, c("gene_name", "DE", "log2FoldChange", "padj")], by = "gene_name")
  
  # Store the DEGs data for the current condition
  degs_in_loops_spec[[cond]] <- degs_uniq
  # Optionally save the DEGs data to file
  #degs_uniq %>% write_tsv(., fs::path(path_results, paste0("/data/", kb, "kb_", toupper(cond), ".DEGs_in_enh_deg_loops.cond_spec.tsv")))
}


# Cluster composition
# Analyze the cluster composition of enhancers involved in condition-specific loops
for(cond in c("scr", "kd")){
  df_loops <- cond_spec[[cond]]
  clusters_in_loops <- c(df_loops$cluster1, df_loops$cluster2)[!is.na(c(df_loops$cluster1, df_loops$cluster2))]
  
  # Print the cluster composition
  print(sprintf("Cluster composition of condition-specific ENH-DEG loops - %s", toupper(cond)))
  print(table(clusters_in_loops))
  cat("\n")
  
  # Calculate high and low cluster counts
  h <- as.data.frame(table(clusters_in_loops)) %>% 
    dplyr::filter(clusters_in_loops %in% clusters_high) %>% 
    dplyr::select(Freq) %>% colSums()
  l <- as.data.frame(table(clusters_in_loops)) %>% 
    dplyr::filter(clusters_in_loops %in% clusters_low) %>% 
    dplyr::select(Freq) %>% colSums()
  
  # Calculate proportions of high and low clusters
  hr <- h / sum(enh_grhl2$group == "high")
  lr <- l / sum(enh_grhl2$group == "low")
  
  # Print proportions and perform chi-square test
  print(sprintf("Number of enh. in clusters HIGH: %s (%s) vs. number of enh. in clusters LOW: %s (%s)", 
                h, round(hr,3), 
                l, round(lr,3)))
  cat("\n")
  
  # Define observed counts for chi-square test
  observed_counts <- matrix(c(h, l, sum(enh_grhl2$group == "high")-h, sum(enh_grhl2$group == "low")-l), nrow = 2, byrow = TRUE)
  colnames(observed_counts) <- c("High", "Low")
  rownames(observed_counts) <- c("In loops", "Not in loops")
  chi_squared_test <- chisq.test(observed_counts)
  
  # Print chi-square test p-value
  print(sprintf("Chi-square test p-value: %s", round(chi_squared_test$p.value,3)))
  cat("\n")
}


##


# GSEA on DEGs-associated loops
library(fgsea)

# Path to Hallmark gene sets for GSEA
hallmark <- path_hallmark

## DEGs in ENH_DEG loops
# Perform GSEA for DEGs in ENH-DEG loops for SCR and KD conditions
gsea_degs <- list()

for(cond in c("scr", "kd")){
  # Rank DEGs by log2FoldChange
  ranked_degs_df <- degs_in_loops[[cond]] %>% arrange(., -log2FoldChange)
  ranked_degs_list <- ranked_degs_df$log2FoldChange
  names(ranked_degs_list) <- ranked_degs_df$gene_name
  
  # Define p-value threshold to filter pathways
  p_val <- 0.05
  
  # Run GSEA
  gsea_degs[[cond]] <- GSEA(gene_list = ranked_degs_list, GO_file = hallmark, p_val = p_val, minSize = 0, maxSize = Inf, collapse = F)
}

# Display GSEA results for SCR and KD
gsea_degs$scr
gsea_degs$kd

# All DEGs (SCR & KD) in ENH-DEG loops
# Combine DEGs from both SCR and KD conditions, rank by log2FoldChange, and run GSEA
ranked_all_df <- rbind(degs_in_loops$scr, degs_in_loops$kd) %>% filter(!duplicated(.)) %>% arrange(., -log2FoldChange)
ranked_all <- ranked_all_df$log2FoldChange
names(ranked_all) <- ranked_all_df$gene_name
gsea_degs_all <- GSEA(gene_list = ranked_all, GO_file = hallmark, p_val = p_val, minSize = 0, maxSize = Inf, collapse = F)
gsea_degs_all


## DEGs in ENH-DEG condition-specific loops
# Perform GSEA for DEGs in ENH-DEG condition-specific loops for SCR and KD conditions
gsea_degs_spec <- list()

for(cond in c("scr", "kd")){
  # Rank DEGs by log2FoldChange
  ranked_degs_df <- degs_in_loops_spec[[cond]] %>% arrange(., -log2FoldChange)
  ranked_degs_list <- ranked_degs_df$log2FoldChange
  names(ranked_degs_list) <- ranked_degs_df$gene_name
  
  # Define p-value threshold to filter pathways
  p_val <- 0.05
  
  # Run GSEA
  gsea_degs_spec[[cond]] <- GSEA(gene_list = ranked_degs_list, GO_file = hallmark, p_val = p_val, minSize = 0, maxSize = Inf, collapse = F)
}

# Display GSEA results for SCR and KD condition-specific loops
gsea_degs_spec$scr
gsea_degs_spec$kd

# All DEGs (SCR & KD) in ENH-DEG condition-specific loops
# Combine DEGs from both SCR and KD condition-specific loops, rank by log2FoldChange, and run GSEA
ranked_all_df <- rbind(degs_in_loops$scr, degs_in_loops_spec$kd) %>% filter(!duplicated(.)) %>% arrange(., -log2FoldChange)
ranked_all <- ranked_all_df$log2FoldChange
names(ranked_all) <- ranked_all_df$gene_name
gsea_degs_all <- GSEA(gene_list = ranked_all, GO_file = hallmark, p_val = p_val, minSize = 0, maxSize = Inf, collapse = F)
gsea_degs_all


##


# Save DEGs for EnrichR
#degs_in_loops_spec[["scr"]] %>% filter(DE == "Down") %>% .$gene_name %>% as.data.frame %>% write_tsv(., fs::path(path_results, "/data/scr_specific_loops_DEGs.Down.tsv"))
#degs_in_loops_spec[["kd"]] %>% filter(DE == "Up") %>% .$gene_name %>% as.data.frame %>% write_tsv(., fs::path(path_results, "/data/kd_specific_loops_DEGs.Up.tsv"))


