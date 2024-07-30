
#
suppressMessages({
  
library(tidyverse)
library(ggplot2)


path_overlaps <- fs::path("/hpcnfs/scratch/P_PGP_FRAGILE_ENHANCERS/results/NEW/data")
path_anno_enhancers <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/results/integrated/annotated_enhancers/NEW/")
path_results_data <- fs::path("/hpcnfs/scratch/P_PGP_FRAGILE_ENHANCERS/results/NEW/data")
path_results_plots <- fs::path("/hpcnfs/scratch/P_PGP_FRAGILE_ENHANCERS/results/NEW/plots")

path_overlaps_icgc <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/results/ICGC/enhancers_SSMs_overlaps/NEW/data/")
#path_overlaps_icgc <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/results/ICGC/enhancers_SSMs_overlaps/data/")
path_results_icgc <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/results/ICGC/enhancers_SSMs_overlaps/NEW/data/")
#path_results_icgc <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/results/ICGC/enhancers_SSMs_overlaps/data/")
  
SEED <- 4321
set.seed(SEED) 
           
WIN <- 3000
MARKERS <- c("CtIP", "GRHL")


##

file_name <- fs::path(path_results_plots, paste0("Rplots.density.snvs_distribution_over_enhancers.pdf")) 
pdf(file_name)

##

# Read enhancers x variants table overlaps 
overlaps <- setNames(vector(mode="list", length=length(MARKERS)), MARKERS)

# FIXME: 
# - overlaps file has duplicated col_names (seqnames, start, end, etc.). 
# - change find_overlaps_snvs.R to save the output with non-duplicated col_names
# - change code below: remove addition of column_names 

for(marker in MARKERS){
  column_names <- c(
    paste0(c("seqnames", "start", "end", "width", "strand", "name", "summit", "cluster"), "_enh"), 
    paste0(c("seqnames", "start", "end", "width", "strand", "REF", "ALT", "QUAL", "FILTER", "SAMPLE", "PURPLE_AF", "AF", "ID"), "_snv")
  )
  overlaps[[marker]] <- read_tsv(fs::path(path_overlaps, paste0(marker, "_enh.hartwig_snvs.overlap.WIN_", WIN, ".tsv")))
  colnames(overlaps[[marker]]) <- column_names
}

# Define groups of clusters
cluster_groups <- setNames(vector(mode="list", length=length(MARKERS)), MARKERS)
cluster_groups$CtIP <- list(
  "high" = c("cluster_1", "cluster_2", "cluster_3", "cluster_4"), 
  "low" = c("cluster_5")
)
cluster_groups$GRHL <- cluster_groups$CtIP # clusters are equal

                     
##


# Plot distribution of SNVs over enhancers regions

anno_only <- F
high_low <- T
each_cluster <- T
save_plots <- F
save_dist <- T
  
for(marker in MARKERS){
  label <- "all_enhancers"
  
  # Distance of SNV from enhancers summit. start_sbj and end_abj are equal
  cat("\n")
  print(paste0("Computing distribution of SNVs distnaces from ", marker, " enhancers summit"))
  
  if(marker == "GRHL" & anno_only == T){
    label <- "annotated_enhancers_only"
    print(paste0("Plotting variants from: ", label))
    annotated_enhancers <- read_tsv(fs::path(path_anno_enhancers, "2kb_GRHL2_enhancers.from_SCR_specific_loops.linked_to_DOWN_DEGs.tsv"))
    overlaps[[marker]] <- overlaps[[marker]] %>% filter(name_enh %in% annotated_enhancers$name)
  }
  
  all_distances <- data.frame("dist" = as.numeric(overlaps[[marker]]$start_snv) - as.numeric(overlaps[[marker]]$summit_enh), 
                              "cluster" = overlaps[[marker]]$cluster_enh, "name" = overlaps[[marker]]$name_enh, 
                              "ID" = overlaps[[marker]]$ID_snv, "PURPLE_AF" = overlaps[[marker]]$PURPLE_AF_snv, "AF" = overlaps[[marker]]$AF_snv)
  if(save_dist == T){
    all_distances %>% write_tsv(., fs::path(path_results_data, paste0("distances.snvs_distribution_over_enhancers.", marker, ".all_clusters.", label, ".tsv")))}
  
  
  # Plot density of SNVs distribution over enhancer region - ALL enhancers
  bw <- 100   
  d_all <- all_distances %>% ggplot(., aes(dist))+
    geom_density(fill = "lightblue3", color = "lightblue3",alpha = 0.7, bw = bw)+
    theme_light()+
    scale_x_continuous(breaks = seq(-3000, 3000, 500))+
    labs(
      title = paste0("SNVs distribution over ", marker, " enhancers"),
      x = "distance from summit", y = ""
    )
  print("Plotting SNVs distribution over enhancer region - ALL enhancers")
  print(d_all)
  
  if(save_plots == T){
  ggsave(plot=d_all, filename=fs::path(path_results_plots, paste0("density.snvs_distribution_over_enhancers.", marker, ".all_clusters.", label, ".png")), 
         device = "png", width = 9.805556, height = 7.027778)}
  

  # Plot density of SNVs distribution over enhancer region - HIGH vs. LOW enhancers
  if(high_low ==T){
    dist_high <- all_distances %>% filter(cluster %in% cluster_groups[[marker]]$high) 
    dist_low <- all_distances %>% filter(cluster %in% cluster_groups[[marker]]$low) 
    
    max_dist_length <- max(dim(dist_high)[1], dim(dist_low)[1])
    df_grouped <- data.frame("high" = c(dist_high$dist, rep(NA, max_dist_length - dim(dist_high)[1])), 
                             "low" = c(dist_low$dist, rep(NA, max_dist_length - dim(dist_low)[1])))
    
    d_high_low <- df_grouped %>% pivot_longer(everything(), names_to = "cluster", values_to = "distances") %>%
      ggplot(., aes(x = distances, color = cluster, fill = cluster))+
      geom_density(bw = bw, alpha = 0.2)+
      theme_light()+
      labs(
        title = paste0("SNVs distribution over ", marker, " enhancers"), 
        x = "distance from summit", y = ""
      )
    print("Plotting SNVs distribution over enhancer region - HIGH vs. LOW enhancers")
    print(d_high_low)
    
    if(save_plots == T){
    ggsave(plot=d_high_low, filename=fs::path(path_results_plots, paste0("density.snvs_distribution_over_enhancers.", marker, ".high_low.", label, ".png")), 
           device = "png", width = 9.805556, height = 7.027778)}
  }
  
  # Plot density of SNVs distribution over enhancer region - EACH cluster of enhancers
  if(each_cluster == T){
    d_each <- all_distances %>% 
      ggplot(., aes(x = dist, gorup = cluster, color = cluster))+
      geom_density(bw = bw, alpha =0.1)+
      theme_light()+
      labs(
        title = paste0("SNVs distribution over ", marker, " enhancers"), 
        x = "distance from summit", y = ""
        )
    print("Plotting SNVs distribution over enhancer region - EACH cluster of enhancers")
    print(d_each)
    
    if(save_plots == T){
    ggsave(plot=d_each, filename=fs::path(path_results_plots, paste0("density.snvs_distribution_over_enhancers.", marker, ".each_cluster.", label, ".png")), 
           device = "png", width = 9.805556, height = 7.027778)}
  }
}


##

dev.off()

##

# TODO: convert into a function 
# TODO: chnage script to accept arguments from command line


##

# Compute & save distances for ICGC 

WIN <- 3000
overlaps_icgc <- setNames(vector(mode="list", length=length(MARKERS)), MARKERS)

for(marker in MARKERS){
  label <- "all_enhancers"
  
  overlaps_icgc[[marker]] <- read_tsv(fs::path(path_overlaps_icgc, paste0("Table_enh_SNVs.", marker, ".all_overlaps.", WIN, "bp_WIN.tsv")))
  overlaps_icgc[[marker]] <- overlaps_icgc[[marker]] %>% dplyr::select(., seqnames, start, end,, width, strand, cluster, summit, name, 
                                                                       seqnames_sbj, start_sbj, end_sbj, width_sbj, strand_sbj, 
                                                                       reference_genome_allele, mutated_to_allele, icgc_donor_id, icgc_sample_id, icgc_mutation_id, AF, ID)
  
  column_names <- c(
    paste0(c("seqnames", "start", "end", "width", "strand", "cluster", "summit", "name"), "_enh"), 
    paste0(c("seqnames", "start", "end", "width", "strand", "REF", "ALT", "DONOR", "SAMPLE", "ICGC_ID", "AF", "ID"), "_snv")
  )
  colnames(overlaps_icgc[[marker]]) <- column_names
  
  all_distances <- data.frame("dist" = as.numeric(overlaps_icgc[[marker]]$start_snv) - as.numeric(overlaps_icgc[[marker]]$summit_enh), 
                              "cluster" = overlaps_icgc[[marker]]$cluster_enh, "name" = overlaps_icgc[[marker]]$name_enh, 
                              "SAMPLE" = overlaps_icgc[[marker]]$SAMPLE_snv, "DONOR" = overlaps_icgc[[marker]]$DONOR_snv, "ID" = overlaps_icgc[[marker]]$ID_snv, 
                              "AF" = overlaps_icgc[[marker]]$AF_snv)
  
  all_distances %>% write_tsv(., fs::path(path_results_icgc, paste0("distances.snvs_distribution_over_enhancers.", marker, ".all_clusters.", label, ".with_AFs.tsv")))
}


##


})
