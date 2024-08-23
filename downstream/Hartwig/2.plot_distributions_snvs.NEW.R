
#
suppressMessages({
  
library(tidyverse)  # Load the tidyverse package for data manipulation and visualization
library(ggplot2)    # Load ggplot2 for creating plots

# Define file paths for input data and results
path_overlaps <- fs::path("/hpcnfs/scratch/P_PGP_FRAGILE_ENHANCERS/results/NEW/data")
path_anno_enhancers <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/results/integrated/annotated_enhancers/NEW/")
path_results_data <- fs::path("/hpcnfs/scratch/P_PGP_FRAGILE_ENHANCERS/results/NEW/data")
path_results_plots <- fs::path("/hpcnfs/scratch/P_PGP_FRAGILE_ENHANCERS/results/NEW/plots")

path_overlaps_icgc <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/results/ICGC/enhancers_SSMs_overlaps/NEW/data/")
path_results_icgc <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/results/ICGC/enhancers_SSMs_overlaps/NEW/data/")

# Set a seed for reproducibility
SEED <- 4321
set.seed(SEED)

# Define window size and markers
WIN <- 3000
MARKERS <- c("CtIP", "GRHL")

##


# Set up the output file for density plots
file_name <- fs::path(path_results_plots, paste0("Rplots.density.snvs_distribution_over_enhancers.pdf")) 
pdf(file_name)


##

# Read and process overlaps between enhancers and SNVs
overlaps <- setNames(vector(mode="list", length=length(MARKERS)), MARKERS)

# FIXME: 
# - overlaps file has duplicated col_names (seqnames, start, end, etc.). 
# - change find_overlaps_snvs.R to save the output with non-duplicated col_names
# - change code below: remove addition of column_names 

for(marker in MARKERS){
  # Define column names for overlaps files
  column_names <- c(
    paste0(c("seqnames", "start", "end", "width", "strand", "name", "summit", "cluster"), "_enh"), 
    paste0(c("seqnames", "start", "end", "width", "strand", "REF", "ALT", "QUAL", "FILTER", "SAMPLE", "PURPLE_AF", "AF", "ID"), "_snv")
  )
  # Read the TSV file for the current marker and assign column names
  overlaps[[marker]] <- read_tsv(fs::path(path_overlaps, paste0(marker, "_enh.hartwig_snvs.overlap.WIN_", WIN, ".tsv")))
  colnames(overlaps[[marker]]) <- column_names
}

# Define cluster groups for each marker
cluster_groups <- setNames(vector(mode="list", length=length(MARKERS)), MARKERS)
cluster_groups$CtIP <- list(
  "high" = c("cluster_1", "cluster_2", "cluster_3", "cluster_4"), 
  "low" = c("cluster_5")
)
cluster_groups$GRHL <- cluster_groups$CtIP  # GRHL has the same cluster groups as CtIP

## 

# Parameters for plotting
anno_only <- FALSE  # Flag to determine if only annotated enhancers should be used
high_low <- TRUE    # Flag to separate plots for high vs. low cluster enhancers
each_cluster <- TRUE  # Flag to create separate plots for each enhancer cluster
save_plots <- FALSE  # Flag to save plots to files
save_dist <- TRUE    # Flag to save distance data to files

# Iterate over each marker to compute and plot SNV distributions
for(marker in MARKERS){
  label <- "all_enhancers"
  
  # Compute distance of SNVs from enhancer summits
  cat("\n")
  print(paste0("Computing distribution of SNVs distances from ", marker, " enhancers summit"))
  
  if(marker == "GRHL" & anno_only){
    label <- "annotated_enhancers_only"
    print(paste0("Plotting variants from: ", label))
    # Read annotated enhancers for GRHL and filter overlaps to include only annotated enhancers
    annotated_enhancers <- read_tsv(fs::path(path_anno_enhancers, "2kb_GRHL2_enhancers.from_SCR_specific_loops.linked_to_DOWN_DEGs.tsv"))
    overlaps[[marker]] <- overlaps[[marker]] %>% filter(name_enh %in% annotated_enhancers$name)
  }
  
  # Create a data frame with distances of SNVs from enhancer summits
  all_distances <- data.frame("dist" = as.numeric(overlaps[[marker]]$start_snv) - as.numeric(overlaps[[marker]]$summit_enh), 
                              "cluster" = overlaps[[marker]]$cluster_enh, "name" = overlaps[[marker]]$name_enh, 
                              "ID" = overlaps[[marker]]$ID_snv, "PURPLE_AF" = overlaps[[marker]]$PURPLE_AF_snv, "AF" = overlaps[[marker]]$AF_snv)
  
  # Save distance data to a file if the flag is set
  if(save_dist){
    all_distances %>% write_tsv(., fs::path(path_results_data, paste0("distances.snvs_distribution_over_enhancers.", marker, ".all_clusters.", label, ".tsv")))
  }
  
  # Plot density of SNVs distribution over enhancer regions for all enhancers
  bw <- 100  # Bandwidth for density plot
  d_all <- all_distances %>% ggplot(., aes(dist)) +
    geom_density(fill = "lightblue3", color = "lightblue3", alpha = 0.7, bw = bw) +
    theme_light() +
    scale_x_continuous(breaks = seq(-3000, 3000, 500)) +
    labs(
      title = paste0("SNVs distribution over ", marker, " enhancers"),
      x = "distance from summit", y = ""
    )
  print("Plotting SNVs distribution over enhancer region - ALL enhancers")
  print(d_all)
  
  # Save the plot if the flag is set
  if(save_plots){
    ggsave(plot=d_all, filename=fs::path(path_results_plots, paste0("density.snvs_distribution_over_enhancers.", marker, ".all_clusters.", label, ".png")), 
           device = "png", width = 9.805556, height = 7.027778)
  }
  
  # Plot density of SNVs distribution for high vs. low enhancer clusters
  if(high_low){
    dist_high <- all_distances %>% filter(cluster %in% cluster_groups[[marker]]$high) 
    dist_low <- all_distances %>% filter(cluster %in% cluster_groups[[marker]]$low) 
    
    # Create a combined data frame for plotting high vs. low clusters
    max_dist_length <- max(dim(dist_high)[1], dim(dist_low)[1])
    df_grouped <- data.frame("high" = c(dist_high$dist, rep(NA, max_dist_length - dim(dist_high)[1])), 
                             "low" = c(dist_low$dist, rep(NA, max_dist_length - dim(dist_low)[1])))
    
    d_high_low <- df_grouped %>% pivot_longer(everything(), names_to = "cluster", values_to = "distances") %>%
      ggplot(., aes(x = distances, color = cluster, fill = cluster)) +
      geom_density(bw = bw, alpha = 0.2) +
      theme_light() +
      labs(
        title = paste0("SNVs distribution over ", marker, " enhancers"), 
        x = "distance from summit", y = ""
      )
    print("Plotting SNVs distribution over enhancer region - HIGH vs. LOW enhancers")
    print(d_high_low)
    
    # Save the plot if the flag is set
    if(save_plots){
      ggsave(plot=d_high_low, filename=fs::path(path_results_plots, paste0("density.snvs_distribution_over_enhancers.", marker, ".high_low.", label, ".png")), 
             device = "png", width = 9.805556, height = 7.027778)
    }
  }
  
  # Plot density of SNVs distribution for each enhancer cluster
  if(each_cluster){
    d_each <- all_distances %>% 
      ggplot(., aes(x = dist, group = cluster, color = cluster)) +
      geom_density(bw = bw, alpha = 0.1) +
      theme_light() +
      labs(
        title = paste0("SNVs distribution over ", marker, " enhancers"), 
        x = "distance from summit", y = ""
      )
    print("Plotting SNVs distribution over enhancer region - EACH cluster of enhancers")
    print(d_each)
    
    # Save the plot if the flag is set
    if(save_plots){
      ggsave(plot=d_each, filename=fs::path(path_results_plots, paste0("density.snvs_distribution_over_enhancers.", marker, ".each_cluster.", label, ".png")), 
             device = "png", width = 9.805556, height = 7.027778)
    }
  }
}

##

# Close the PDF device
dev.off()

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
