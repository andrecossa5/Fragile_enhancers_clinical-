
suppressMessages({
  

library(tidyverse)
library(ggplot2)
library(ggpubr)

MARKERS <- c("CtIP", "GRHL")
SOURCES <- c("hart", "icgc")
anno_only <- F
label <- "all_enhancers"

path_anno_enhancers <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/results/integrated/annotated_enhancers/NEW/")
hart.ctip.dist <- read_tsv("/hpcnfs/scratch/P_PGP_FRAGILE_ENHANCERS/results/NEW/data/distances.snvs_distribution_over_enhancers.CtIP.all_clusters.all_enhancers.tsv")
hart.grhl.dist <- read_tsv("/hpcnfs/scratch/P_PGP_FRAGILE_ENHANCERS/results/NEW/data/distances.snvs_distribution_over_enhancers.GRHL.all_clusters.all_enhancers.tsv")
icgc.ctip.dist <- read_tsv("/hpcnfs/scratch/PGP/Ciacci_et_al/results/ICGC/enhancers_SSMs_overlaps/NEW/data/distances.snvs_distribution_over_enhancers.CtIP.all_clusters.all_enhancers.with_AFs.tsv")
icgc.grhl.dist <- read_tsv("/hpcnfs/scratch/PGP/Ciacci_et_al/results/ICGC/enhancers_SSMs_overlaps/NEW/data/distances.snvs_distribution_over_enhancers.GRHL.all_clusters.all_enhancers.with_AFs.tsv")
annotated_enhancers <- read_tsv(fs::path(path_anno_enhancers, "2kb_GRHL2_enhancers.from_SCR_specific_loops.linked_to_DOWN_DEGs.tsv")) %>% suppressMessages()

path_results_plots <- fs::path("/hpcnfs/scratch/P_PGP_FRAGILE_ENHANCERS/results/NEW/plots/")
path_results_data <- fs::path("/hpcnfs/scratch/P_PGP_FRAGILE_ENHANCERS/results/NEW/data/")

# Name of pdf files with plots
file_name1 <- fs::path(path_results_plots, paste0("Rplots.3.distributions_across_clusters.pdf")) 
file_name2 <- fs::path(path_results_plots, paste0("Rplots.3.distributions_across_clusters.part_2.pdf")) 

if(anno_only == T){
  label <- "only_annotated_ehancers_GRHL2"
  hart.grhl.dist <- hart.grhl.dist %>% filter(name %in% annotated_enhancers$name)
  icgc.grhl.dist <- icgc.grhl.dist %>% filter(name %in% annotated_enhancers$name)
  file_name1 <- fs::path(path_results_plots, paste0("Rplots.3.distributions_across_clusters.anno_only.pdf")) 
  file_name2 <- fs::path(path_results_plots, paste0("Rplots.3.distributions_across_clusters.part_2.anno_only.pdf")) 
}

all_dfs <- list("CtIP" = list("hart" = hart.ctip.dist, "icgc" = icgc.ctip.dist), 
                "GRHL" = list("hart" = hart.grhl.dist, "icgc" = icgc.grhl.dist))


##

pdf(file_name1, width = 11, height = 7)

## 


max_length <- max(sapply(list(hart.ctip.dist, hart.grhl.dist, icgc.ctip.dist, icgc.grhl.dist), function(x){dim(x)[1]}))

# Densities hartwig vs. icgc - ALL enhancers
df_full <- data.frame(
  "hart.ctip" = c(hart.ctip.dist$dist, rep(NA, max_length-dim(hart.ctip.dist)[1])), 
  "hart.grhl" = c(hart.grhl.dist$dist, rep(NA, max_length-dim(hart.grhl.dist)[1])),
  "icgc.ctip" = c(icgc.ctip.dist$dist, rep(NA, max_length-dim(icgc.ctip.dist)[1])),
  "icgc.grhl" = c(icgc.grhl.dist$dist, rep(NA, max_length-dim(icgc.grhl.dist)[1]))
)


##


for(marker in tolower(MARKERS)){
  df <- df_full[, endsWith(colnames(df_full), marker)]
  d <- df %>% pivot_longer(everything(), names_to = "source", values_to = "distances") %>%
    ggplot(., aes(x = distances, color = source, fill = source))+
    geom_density(bw = 100, alpha = 0.2)+
    theme_light()+
    labs(title = paste0("All clusters of enhancers - ", marker), 
                        subtitle = label)
  print(d)
}


##


# Define groups of clusters
cluster_groups <- setNames(vector(mode="list", length=length(MARKERS)), MARKERS)
cluster_groups$CtIP <- list(
  "high" = c("cluster_1", "cluster_2", "cluster_3", "cluster_4"), 
  "low" = c("cluster_5")
)
cluster_groups$GRHL <- cluster_groups$CtIP # clusters are equal


# x cluster high-low
df_high <- data.frame(
  "hart.ctip" = c(hart.ctip.dist[hart.ctip.dist$cluster %in% cluster_groups$CtIP$high, ]$dist, rep(NA, max_length-dim(hart.ctip.dist[hart.ctip.dist$cluster %in% cluster_groups$CtIP$high, ])[1])), 
  "hart.grhl" = c(hart.grhl.dist[hart.grhl.dist$cluster %in% cluster_groups$GRHL$high, ]$dist, rep(NA, max_length-dim(hart.grhl.dist[hart.grhl.dist$cluster %in% cluster_groups$GRHL$high, ])[1])),
  "icgc.ctip" = c(icgc.ctip.dist[icgc.ctip.dist$cluster %in% cluster_groups$CtIP$high, ]$dist, rep(NA, max_length-dim(icgc.ctip.dist[icgc.ctip.dist$cluster %in% cluster_groups$CtIP$high, ])[1])),
  "icgc.grhl" = c(icgc.grhl.dist[icgc.grhl.dist$cluster %in% cluster_groups$GRHL$high, ]$dist, rep(NA, max_length-dim(icgc.grhl.dist[icgc.grhl.dist$cluster %in% cluster_groups$GRHL$high, ])[1]))
)
df_low <- data.frame(
  "hart.ctip" = c(hart.ctip.dist[hart.ctip.dist$cluster %in% cluster_groups$CtIP$low, ]$dist, rep(NA, max_length-dim(hart.ctip.dist[hart.ctip.dist$cluster %in% cluster_groups$CtIP$low, ])[1])), 
  "hart.grhl" = c(hart.grhl.dist[hart.grhl.dist$cluster %in% cluster_groups$GRHL$low, ]$dist, rep(NA, max_length-dim(hart.grhl.dist[hart.grhl.dist$cluster %in% cluster_groups$GRHL$low, ])[1])),
  "icgc.ctip" = c(icgc.ctip.dist[icgc.ctip.dist$cluster %in% cluster_groups$CtIP$low, ]$dist, rep(NA, max_length-dim(icgc.ctip.dist[icgc.ctip.dist$cluster %in% cluster_groups$CtIP$low, ])[1])),
  "icgc.grhl" = c(icgc.grhl.dist[icgc.grhl.dist$cluster %in% cluster_groups$GRHL$low, ]$dist, rep(NA, max_length-dim(icgc.grhl.dist[icgc.grhl.dist$cluster %in% cluster_groups$GRHL$low, ])[1]))
)

# High 
for(marker in tolower(MARKERS)){
  df <- df_high[, endsWith(colnames(df_high), marker)]
  d <- df %>% pivot_longer(everything(), names_to = "source", values_to = "distances") %>%
    ggplot(., aes(x = distances, color = source, fill = source))+
    geom_density(bw = 100, alpha = 0.2)+
    theme_light()+
    labs(title = paste0("High clusters - ", marker), 
         subtitle = label)
  print(d)
}


# Low
for(marker in tolower(MARKERS)){
  df <- df_low[, endsWith(colnames(df_low), marker)]
  d <- df %>% pivot_longer(everything(), names_to = "source", values_to = "distances") %>%
    ggplot(., aes(x = distances, color = source, fill = source))+
    geom_density(bw = 100, alpha = 0.2)+
    theme_light()+
    labs(title = paste0("Low clusters - ", marker), 
         subtitle = label)
  print(d)
}

##


# Each cluster 

plot_list <- list(
  "CtIP" = setNames(vector(mode="list", length=length(unlist(cluster_groups[["CtIP"]]))), unlist(cluster_groups[["CtIP"]])),
  "GRHL" = setNames(vector(mode="list", length=length(unlist(cluster_groups[["GRHL"]]))), unlist(cluster_groups[["GRHL"]]))
)

for(marker in MARKERS){
  clusts <- unlist(cluster_groups[[marker]])
  for(clust in clusts){
    print(clust)
    df1 <- all_dfs[[marker]][["hart"]] %>% filter(cluster == clust)
    df2 <- all_dfs[[marker]][["icgc"]] %>% filter(cluster == clust)
    
    max_length <- max(dim(df1)[1], dim(df2)[1])
    df <- data.frame("hart" = c(df1$dist, rep(NA, max_length-length(df1$dist))), 
                     "icgc" = c(df2$dist, rep(NA, max_length-length(df2$dist))))
    
    d <- df %>% pivot_longer(everything(), names_to = "source", values_to = "distances") %>%
      ggplot(., aes(x=distances, group = source, color = source))+
      geom_density(bw=100, alpha=0.1)+
      theme_light()+
      labs(
        title = paste0(clust, " ", marker), 
        subtitle = label, 
        x = "distance from summit", y = "")
    print(d)
    
    plot_list[[marker]][[clust]] <- d
    }
}

#for(marker in (MARKERS)){
  #cowplot::plot_grid("title", plotlist = plot_list[[marker]], ncol = 1) #%>%
    #ggsave(., filename = paste0("/Users/ieo6983/Desktop/fragile_enhancer_clinical/results/Hartwig/snvs_distribution_over_enhancers/plots/density.snvs_distribution_over_enhancers.", marker,".each_cluster.hartwig_vs_icgc.", label, ".png"), device = "png", width = 7, height = 22)}


##


### ANALYSES on VAFs ###

WIN <- 100

df_dist_wins <- list(
  "CtIP" = list("hart", "icgc"), 
  "GRHL" = list("hart", "icgc")
)
  
snvs_x_enh <- list(
  "CtIP" = list("hart", "icgc"), 
  "GRHL" = list("hart", "icgc")
)

df_n_muts_occ <- list(
  "CtIP" = list("hart", "icgc"), 
  "GRHL" = list("hart", "icgc")
)


for(source in SOURCES){
  cat("\n")
  print(paste0(" ----- ", source, " ----- "))
  print("Storing variants information")
  for(marker in MARKERS){
    df.dist.win <- all_dfs[[marker]][[source]] %>% filter(abs(dist) <= WIN)
    df_dist_wins[[marker]][[source]] <- df.dist.win
    
    # Store number of SNVs x enhancer info
    n_muts_x_enh <- df.dist.win %>% group_by(., name) %>% summarise("n_muts_x_enh" = n())
    n_muts_x_enh <- left_join(n_muts_x_enh, df.dist.win[, c("cluster", "name")], by = "name")
    snvs_x_enh[[marker]][[source]] <- n_muts_x_enh
    
    # Store number of occurrences x SNV info
    n_muts_occ <- df.dist.win %>% group_by(., ID) %>% summarise(n_muts_occ = n())
    df_n_muts_occ[[marker]][[source]] <- n_muts_occ
  }
}


# PLOTS
for(source in SOURCES){
  cat("\n")
  print(paste0(" ----- ", source, " ----- "))
  for(marker in MARKERS){
    df.dist.win <- df_dist_wins[[marker]][[source]]
    
    # Plot number of mutations x enhancer across donors
    sort(table(df.dist.win$name), decreasing = T) %>% 
      barplot(las=3, cex.names = 0.5, 
              main=paste0("Number of SNVs x ", marker ," enhancer across all ", source, " samples: ", WIN, " bp win"))
    up_4 <- table(df.dist.win$name)[table(df.dist.win$name) > 4]
    up_6 <- table(df.dist.win$name)[table(df.dist.win$name) > 6]
    print("--- Number of SNVs per enhancer ---")
    print(paste0("Number of ", marker, " enhancers with > 4 SNVs: ", length(up_4)))
    print(paste0("Number of ", marker, " enhancers with > 6 SNVs: ", length(up_6)))
    print(up_6)
    cat("\n")
    
    # Print number of occurrences of each SNV across donors
    sort(table(df.dist.win$ID), decreasing = T) %>% 
      barplot(las=3, cex.names = 0.5,  
              main=paste0("Number of occurrences of SNVs in ", marker ," enhancer across all ", source, " samples: ", WIN, " bp win"))
    up_1 <- table(df.dist.win$ID)[table(df.dist.win$ID) > 1]
    up_2 <- table(df.dist.win$ID)[table(df.dist.win$ID) > 2]
    print("--- Number of recurring variants ---")
    print(paste0("Number of SNVs with > 1 occ. in ", marker, " enhancers: ", length(up_1)))
    print(paste0("Number of SNVs with > 2 occ. in ", marker, " enhancers: ", length(up_2)))
    
    # Merge info
    n_muts_x_enh <- snvs_x_enh[[marker]][[source]]
    n_muts_occ <- df_n_muts_occ[[marker]][[source]]
    all_occurrences <- left_join(df.dist.win, n_muts_x_enh, by = "name") %>% 
      left_join(., n_muts_occ, by = "ID")
    
    # Number of enhancers with SNVs with occurrence > ?
    n_enh <- length(unique(all_occurrences[all_occurrences$n_muts_x_enh > 2, ]$name))
    print(paste0("High-occurrence SNVs in ", source, " corresponding to ", n_enh, " ", marker, " enhancers"))
    cat("\n")
    
    # Plot AF distribution of SNVs in enhancers divided by number of occurrences    
    af <- all_occurrences %>%
      ggplot(., aes(x = n_muts_occ, y = AF))+
      geom_point()+
      theme_light()+
      labs(x = "N. of SNV occurrences", 
           title = paste0("AF of high-occurrence SNVs in ", source, " - ", marker))
    print(af)
    
    aff <- all_occurrences %>% filter(., n_muts_occ > 2) %>% 
      ggplot(., aes(x = n_muts_occ, y = AF))+
      geom_point()+
      theme_light()+
      labs(x = "N. of SNV occurrences", 
           title = paste0("AF of high-occurrence SNVs in ", source, " - ", marker))
    print(aff)
  }
}


##

common_enh_sub_dfs <- list("CtIP", "GRHL")
for(marker in MARKERS){
  df.dist.marker <- df_dist_wins[[marker]]
  
  AF_i <- df.dist.marker[["icgc"]][, c("name", "cluster", "AF")]
  AF_i$source <- rep("icgc", dim(AF_i)[1])
  AF_i$group <- "high"
  AF_i[AF_i$cluster %in% cluster_groups[[marker]]$low, ]$group <- "low"
  
  AF_h <- df.dist.marker[["hart"]][, c("name", "cluster", "AF")]
  AF_h$source <- rep("hart", dim(AF_h)[1])
  AF_h$group <- "high"
  AF_h[AF_h$cluster %in% cluster_groups[[marker]]$low, ]$group <- "low"
  
  AF_all <- na.omit(rbind(AF_i, AF_h)) # Many ICGC samples have no AF info 
  AF_all$source <- factor(AF_all$source, levels = c("icgc", "hart"))

  # Plot AFs distributions for all enhancers 
  p_all <- AF_all %>% 
    ggplot(., aes(x = source, y = AF, fill = source))+
    geom_boxplot()+
    stat_compare_means(label = "p.signif", method = "wilcox.test", 
                       size = 3, comparisons = list(c("icgc", "hart")))+
    scale_fill_manual(values = c("icgc" = "#9fc8c8", "hart" = "#298c8c"))+
    theme_light()+
    labs(title = paste0("AFs distribution of all SNVs within enhancers ", marker, " - icgc vs. hartwig"), 
         subtitle = "all enhancers")
  print(p_all)
  
  p_groups <- AF_all %>% 
    ggplot(., aes(x = source, y = AF, fill = source))+
    geom_boxplot()+
    facet_wrap(~group)+
    stat_compare_means(label = "p.signif", method = "wilcox.test", 
                       size = 3, comparisons = list(c("icgc", "hart")))+
    scale_fill_manual(values = c("icgc" = "#9fc8c8", "hart" = "#298c8c"))+
    theme_light()+
    labs(title = paste0("AFs distribution of all SNVs within enhancers ", marker, " - icgc vs. hartwig"), 
         subtitle = "high vs. low")
  print(p_groups)
  
  p_clusters <- AF_all %>% 
    ggplot(., aes(x = source, y = AF, fill = source))+
    geom_boxplot()+
    facet_wrap(~cluster, nrow = 1)+
    stat_compare_means(label = "p.signif", method = "wilcox.test", 
                       size = 3, comparisons = list(c("icgc", "hart")))+
    scale_fill_manual(values = c("icgc" = "#9fc8c8", "hart" = "#298c8c"))+
    theme_light()+
    labs(title = paste0("AFs distribution of all SNVs within enhancers ", marker, " - icgc vs. hartwig"), 
         subtitle = "x cluster")
  print(p_clusters)
  
  ###
  
  # Compute average AF of SNVs for each enhancer 
  avg_AF_x_enh <- AF_all %>% group_by(source, cluster, name) %>% dplyr::summarise(avg_af = mean(AF))
  splitted_dfs_list <- split(avg_AF_x_enh, avg_AF_x_enh$source)
  merged_avg_AF <- full_join(splitted_dfs_list[[1]], splitted_dfs_list[[2]], by = "name")
  
  tot_enh <- length(unique(merged_avg_AF$name)) 
  print(paste0("Total number of enhancers with SNVs in icgc and hart (union):", tot_enh))
  
  common_enh <- merged_avg_AF %>% filter(., !is.na(source.x) & ! is.na(source.y))
  print(paste0("Number of enhancers with SNVs in icgc and hart (both):", length(unique(common_enh$name))))
  perc <- round(length(unique(common_enh$name)) / tot_enh * 100, 2)
  print( paste0( perc , " %"))
    
  # Compute FC of avg_af in hartwig vs. icgc, for common enhancers 
  common_enh$avg_af.fc <- log2(common_enh$avg_af.y / common_enh$avg_af.x)
  t_test_res <- t.test(common_enh$avg_af.fc)$p.value
  fc_hist <- common_enh %>% ggplot(., aes(x = avg_af.fc))+
    geom_histogram()+
    theme_light()+
    scale_x_continuous(n.breaks = 10)+
    geom_vline(xintercept = mean(common_enh$avg_af.fc), color = "red")+
    geom_vline(xintercept = median(common_enh$avg_af.fc), color = "red4")+
    labs(x = "Avg. AF log2(Fold-change)", y = "", 
         title = "Distribution of avg.AF fold-change in hartwig vs. icgc", 
         subtitle = paste0("For ", length(unique(common_enh$name)), " ", marker, " common enhancers (", perc, "%)"), 
         caption = paste0("T-test for mean different from 0 - P-value =  ", round(t_test_res, 3)))
  print(fc_hist)
  
  # Common enhancers identity 
  # High - Low
  common_enh$group <- "high"
  common_enh[common_enh$cluster.x %in% cluster_groups[[marker]]$low, ]$group <- "low"
  plot_pie1 <- common_enh %>% group_by(., group) %>% summarise(., counts = n())
  pie_groups <- ggplot(plot_pie1, aes(x="", y=counts, fill=group)) +
    geom_bar(stat="identity", width=1, color = "white") +
    coord_polar("y", start=0) +
    theme_void()
  print(pie_groups)
  
  # Each cluster
  plot_pie2 <- common_enh %>% group_by(., cluster.x) %>% summarise(., counts = n())
  pie_clusters <- ggplot(plot_pie2, aes(x="", y=counts, fill=cluster.x)) +
    geom_bar(stat="identity", width=1, color = "white") +
    coord_polar("y", start=0) +
    theme_void()
  print(pie_clusters)
  
  # Annotated (loops) or not
  common_enh$anno <- "no"
  common_enh[common_enh$name %in% annotated_enhancers$name, ]$anno <- "yes"
  print("Which annotated enhancers: ")
  print(common_enh[common_enh$anno == "yes", ]$name)
  plot_pie3 <- common_enh %>% group_by(., anno) %>% summarise(., counts = n())
  pie_anno <- ggplot(plot_pie3, aes(x="", y=counts, fill=anno)) +
    geom_bar(stat="identity", width=1, color = "white") +
    coord_polar("y", start=0) +
    theme_void()
  print(pie_anno)
  # Save df of annotated enhancers with mutations across datasets 
  common_enh[common_enh$anno == "yes", ] %>% write_tsv(., 
            fs::path(path_results_data, paste0("Anno_enhancers.mutated_across_hart_and_icgc.", marker, ".tsv")))
  
  # Only common enhancers wit log2FC(avg. AF) > 0 
  common_enh_sub <- common_enh %>% filter(., avg_af.fc > 0)
  common_enh_sub_dfs[[marker]] <- common_enh_sub
  plot_pie1 <- common_enh_sub %>% group_by(., group) %>% summarise(., counts = n())
  pie_groups <- ggplot(plot_pie1, aes(x="", y=counts, fill=group)) +
    geom_bar(stat="identity", width=1, color = "white") +
    coord_polar("y", start=0) +
    theme_void()
  print(pie_groups)
  plot_pie2 <- common_enh_sub %>% group_by(., cluster.x) %>% summarise(., counts = n())
  pie_clusters <- ggplot(plot_pie2, aes(x="", y=counts, fill=cluster.x)) +
    geom_bar(stat="identity", width=1, color = "white") +
    coord_polar("y", start=0) +
    theme_void()
  print(pie_clusters)
  
}


##

dev.off()

##


#

pdf(file_name2, width = 15, height = 7)

#

### SNVs conservation across datasets and VAFs 

for(marker in MARKERS){
  df.dist.marker <- df_dist_wins[[marker]]
  
  AF_i <- df.dist.marker[["icgc"]][, c("name", "cluster", "ID", "AF")]
  AF_i$source <- rep("icgc", dim(AF_i)[1])
  AF_i$group <- "high"
  AF_i[AF_i$cluster %in% cluster_groups[[marker]]$low, ]$group <- "low"
  
  AF_h <- df.dist.marker[["hart"]][, c("name", "cluster", "ID", "AF")]
  AF_h$source <- rep("hart", dim(AF_h)[1])
  AF_h$group <- "high"
  AF_h[AF_h$cluster %in% cluster_groups[[marker]]$low, ]$group <- "low"
  
  AF_all <- na.omit(rbind(AF_i, AF_h)) # Many ICGC samples have no AF info 
  AF_all$source <- factor(AF_all$source, levels = c("icgc", "hart"))
  
  
  # SNVs in common
  AF_all$POS_ID <- str_sub(str_sub(AF_all$ID, end = -5))
  splitted_dfs_list <- split(AF_all, AF_all$source)
  
  ## Level 1 - SNVs equal for position 
  tot_SNVs_pos <- length(unique(AF_all$POS_ID)) 
  print(paste0("Total number of POS-SNVs in icgc and hart (union):", tot_SNVs_pos))
  common_SNVs_pos <- splitted_dfs_list$icgc$POS_ID[splitted_dfs_list$icgc$POS_ID %in% splitted_dfs_list$hart$POS_ID]
  common_SNVs_pos_df <- AF_all %>% dplyr::filter(., POS_ID %in% common_SNVs_pos)
  print(paste0("Number of SNVs in icgc and hart (both):", length(unique(common_SNVs_pos_df$POS_ID))))
  
  if(length(unique(common_SNVs_pos_df$POS_ID)) != 0){
    # Occurrences Level 1
    df_occ_pos <- common_SNVs_pos_df %>% group_by(POS_ID, source) %>% dplyr::summarise(., n_occ = n()) %>% ungroup()
    df_occ_pos$source <- factor(df_occ_pos$source, levels = c("icgc", "hart")) 
    p_occ1 <- df_occ_pos %>%
      ggplot(., aes(x=source, y = n_occ, fill = source))+
      geom_bar(stat = "identity", position = "dodge")+
      facet_wrap(~POS_ID, nrow = 1)+
      labs(x = "", y = "number of occurrences", 
           title = paste0("Number of occurrences for common SNVs positions - ", marker))+
      theme(strip.text = element_text(size = 7))
    print(p_occ1)
    
    # AFs of SNVs in common - Level 1 - POS
    common_SNVs_pos_df$AF <- as.numeric(common_SNVs_pos_df$AF) 
    common_SNVs_pos_df$source <- factor(common_SNVs_pos_df$source, levels = c("icgc", "hart"))
    
    p_af1 <- common_SNVs_pos_df %>%
      ggplot(., aes(x = source, y = AF, fill = source))+
      geom_segment(aes(x=source, xend=source, y=0, yend=AF), color="grey")+
      geom_point( color="orange", size=4)+
      facet_wrap(~POS_ID, nrow = 1)+
      theme_light()+
      labs(x = "",
           title = paste0("AF values for common POS-SNVs - ", marker))+
      theme(strip.text = element_text(size = 7))
    print(p_af1)
    
    ## Level 2 - SNVs equal for position AND base substitution
    tot_SNVs <- length(unique(AF_all$ID)) 
    print(paste0("Total number of SNVs in icgc and hart (union):", tot_SNVs))
    common_SNVs <- splitted_dfs_list$icgc$ID[splitted_dfs_list$icgc$ID %in% splitted_dfs_list$hart$ID]
    common_SNVs_df <- AF_all %>% dplyr::filter(., ID %in% common_SNVs)
    print(paste0("Number of SNVs in icgc and hart (both):", length(unique(common_SNVs_df$ID))))
    
    # Occurrences Level 2
    df_occ <- common_SNVs_df %>% group_by(ID, source) %>% dplyr::summarise(., n_occ = n()) %>% ungroup()
    df_occ$source <- factor(df_occ$source, levels = c("icgc", "hart")) 
    p_occ2 <- df_occ %>%
      ggplot(., aes(x=source, y = n_occ, fill = source))+
      geom_bar(stat = "identity", position = "dodge")+
      facet_wrap(~ID, nrow = 1)+
      labs(x = "", y = "number of occurrences", 
           title = paste0("Number of occurrences for common SNVs - ", marker))+
      theme(strip.text = element_text(size = 7))
    print(p_occ2)
    
    ## AFs of SNVs in common - Level 2 - POS-BASE
    common_SNVs_df$AF <- as.numeric(common_SNVs_df$AF) 
    common_SNVs_df$source <- factor(common_SNVs_df$source, levels = c("icgc", "hart"))
    
    p_af2 <- common_SNVs_df %>%
      ggplot(., aes(x = source, y = AF, group = source))+
      geom_segment(aes(x=source, xend=source, y=0, yend=AF), color="grey")+
      geom_point( color="orange", size=4)+
      facet_wrap(~POS_ID, nrow = 1)+
      theme_light()+
      labs(x = "",
           title = paste0("AF values for common SNVs - ", marker))+
      theme(strip.text = element_text(size = 7))
    print(p_af2)
  }
  
  ## All SNVs in common_enhancers with log2FC(avg. AF) > 0
  AFs_common_enh <- AF_all %>% dplyr::filter(name %in% common_enh_sub_dfs[[marker]]$name) 
  p_hist <- AFs_common_enh %>% filter(., source == "hart") %>%
    ggplot(., aes(x = AF))+
    geom_histogram()+
    theme_light()+
    scale_x_continuous(n.breaks = 10)+
    labs(x = "AF", y = "", 
         title = paste0("Distribution of hartwig AFs for SNVs within enhancers - ", marker), 
         subtitle = "Only enhancers with log2FC(avg. AF) > 0")
  print(p_hist)
  
  # Plot in detail only top SNVs in terms of AF
  p_af_top_common <- AFs_common_enh %>% filter(., source == "hart") %>% 
    arrange(., -AF) %>% top_n(n=25,wt=AF) %>%
    ggplot(., aes(x = ID, y = AF, group = ID))+
    geom_segment(aes(x=ID, xend=ID, y=0, yend=AF), color = "grey")+
    geom_point(color = "#a00000", size = 4)+
    facet_wrap(~name, nrow=1, scales = "free_x")+
    theme_light()+
    labs(title = paste0("Hartwig AF values for SNVs within enahcners - ", marker), 
         subtitle = "Only enhancers with log2FC(avg. AF) > 0")+
    theme(strip.text = element_text(size = 5), 
          axis.text.x = element_text(angle = 90, vjust = 0.5))
  print(p_af_top_common)
  
  # Save table
  table_common <- AFs_common_enh %>% filter(., source == "hart") %>% 
    arrange(., -AF)
  table_common %>% write_tsv(., fs::path(path_results_data, paste0("SNVs_top_AFs.log2FC(avgAF)_greater_0.common_enhancers_", marker, ".tsv")))
  
  
  ##
  
  # Types of substitutions
  table_common$sub <- str_sub(table_common$ID, start = -3)
  SBS_types <- table_common %>% group_by(., sub) %>%
    dplyr::summarise(n = n()) %>%
    ggplot(., aes(x = sub, y = n, fill = sub))+
    geom_bar(stat = "identity", position = "dodge")+
    theme_light()+
    labs(x = "", y = "", 
          title = paste0("Types of SBSs - ", marker))
  print(SBS_types)
  
  # Save variants within common enahncers for mut.sig.analysis
  common_enh_var_df <- AF_all[AF_all$name %in% common_enh_sub_dfs[[marker]]$name,]
  common_enh_var_df %>% write_tsv(., fs::path(path_results_data, paste0("SNVs_coords.common_enhancers_mutated_across_hart_and_icgc.", marker, ".tsv")))
  
}


##

dev.off()

##

})

##



