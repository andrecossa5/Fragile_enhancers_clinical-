
library(tidyverse)
library(ggpubr)
library(rowr)

SEED <- 4321
set.seed(SEED)

path_output_temp <- fs::path("/hpcnfs/scratch/P_PGP_FRAGILE_ENHANCERS/results/temp/")

SOURCES <- c("icgc", "hart")
MARKERS <- c("CtIP", "GRHL")
WIN <- 50
Mb <- 1000000

high_low <- T


##


# Define groups of clusters
cluster_groups <- setNames(vector(mode="list", length=length(MARKERS)), MARKERS)
cluster_groups$CtIP <- list(
  "high" = c("CtIP_cluster_2_Enh", "CtIP_cluster_3_Enh", "CtIP_cluster_5_Enh"), 
  "low" = c("CtIP_cluster_6.2_Enh", "CtIP_cluster_6.3_Enh", "CtIP_cluster_6.1_Enh", "CtIP_cluster_6.0")
)
cluster_groups$GRHL <- list(
  "high" = c("GRHL_cluster_1_Enh", "GRHL_cluster_2_Enh", "GRHL_cluster_4_Enh"), 
  "low" = c("GRHL_cluster_5.2_Enh", "GRHL_cluster_5.3_Enh", "GRHL_cluster_5.1_Enh", "GRHL_cluster_5.0")
)

# Read FC distributions from ICGC and Hartwig
fc_all_samples_markers_source <- list()
fc_ran_all_samples_markers_source <- list()
fc_x_group_all_samples_markers_source <- list()
for(source in SOURCES){
  fc_all_samples_markers_source[[source]] <- readRDS(file = fs::path(path_output_temp, paste0("fc_dist.all_enhancers.", source, ".rds")))   
  fc_ran_all_samples_markers_source[[source]] <- readRDS(file = fs::path(path_output_temp, paste0("fc_dist.all_random.", source, ".rds")))  
  fc_x_group_all_samples_markers_source[[source]] <- readRDS(file = fs::path(path_output_temp, paste0("fc_dist.grouped_enhancers.", source, ".rds")))   
}


##

pdf_file_name <- fs::path("/hpcnfs/scratch/P_PGP_FRAGILE_ENHANCERS/results/temp/Rplots.mutation_rate_icgc_vs_hartwig.all_enhancers.pdf")
pdf(pdf_file_name)

##


## PLOT

for(marker in MARKERS){
  
  # Mutation freq. (FC over backround mut. rate) in ICGC vs. Hartwig, for all enhancers
  
  df_to_plot <- cbind.fill(fc_all_samples_markers_source[["icgc"]][[marker]], 
             fc_all_samples_markers_source[["hart"]][[marker]], 
             fc_ran_all_samples_markers_source[["icgc"]][[marker]],
             fc_ran_all_samples_markers_source[["hart"]][[marker]])
  colnames(df_to_plot) <- c("fc_enh_icgc", "fc_enh_hart", "fc_ran_icgc", "fc_ran_hart")
    
  df_to_plot_long <- df_to_plot %>% pivot_longer(everything(), names_to = "seq_type", values_to = "mut_rate")
  df_to_plot_long$seq_type <- factor(df_to_plot_long$seq_type, levels = c("fc_enh_icgc", "fc_enh_hart", "fc_ran_icgc", "fc_ran_hart"))
  p <- df_to_plot_long %>% 
    ggplot(., aes(x = seq_type, y = mut_rate, group = seq_type))+
    geom_jitter(width = 0.3, alpha = 0.4)+
    geom_boxplot(alpha=0.2, width = 0.4, outliers = F)+
    stat_compare_means(label = "p.signif", method = "wilcox.test", size = 3, 
                       comparisons = list(c("fc_enh_icgc", "fc_enh_hart")))+
    theme_light()+
    geom_hline(yintercept = 1, col = "grey", linetype = 2)+
    labs(x = "", y = "Fold-change over background", 
         title = paste0("Mutation frequency in ", marker, " enhancers"), 
         subtitle = "ICGC vs. Hartwig")+
    stat_summary(
      fun = mean, 
      geom = "text", 
      aes(label = sprintf("%.2f", ..y..)), 
      vjust = -20,
      color = "blue"
    )
  print(p)
  
  # Same jitter_plot, but in log2 scale
  
  df_to_plot_log <- cbind.fill(log2(fc_all_samples_markers_source[["icgc"]][[marker]]+0.1), 
                           log2(fc_all_samples_markers_source[["hart"]][[marker]]+0.1), 
                           log2(fc_ran_all_samples_markers_source[["icgc"]][[marker]]+0.1),
                           log2(fc_ran_all_samples_markers_source[["hart"]][[marker]]+0.1))
  colnames(df_to_plot_log) <- c("fc_enh_icgc", "fc_enh_hart", "fc_ran_icgc", "fc_ran_hart")
  df_to_plot_log_long <- df_to_plot_log %>% pivot_longer(everything(), names_to = "seq_type", values_to = "mut_rate")
  df_to_plot_log_long$seq_type <- factor(df_to_plot_log_long$seq_type, levels = c("fc_enh_icgc", "fc_enh_hart", "fc_ran_icgc", "fc_ran_hart"))
  
  p <- df_to_plot_log_long %>% 
    ggplot(., aes(x = seq_type, y = mut_rate, group = seq_type))+
    geom_jitter(width = 0.3, alpha = 0.4)+
    geom_boxplot(alpha=0.2, width = 0.4, outliers = F)+
    stat_compare_means(label = "p.signif", method = "wilcox.test", size = 3, 
                       comparisons = list(c("fc_enh_icgc", "fc_enh_hart")))+
    theme_light()+
    geom_hline(yintercept = 0, col = "grey", linetype = 2)+
    labs(x = "", y = "log2(Fold-change over background)", 
         title = paste0("Mutation frequency in ", marker, " enhancers"), 
         subtitle = "ICGC vs. Hartwig")+
    stat_summary(
      fun = mean, 
      geom = "text", 
      aes(label = sprintf("%.2f", ..y..)), 
      vjust = -20,
      color = "blue"
    )
  print(p)
  
  
  ##
  
  
  # Per group - high vs. low
  
  df_to_plot <- cbind.fill(fc_all_samples_markers_source[["icgc"]][[marker]], 
                           fc_all_samples_markers_source[["hart"]][[marker]],
                           fc_x_group_all_samples_markers_source[["icgc"]][[marker]][["high"]], # high ICGC 
                           fc_x_group_all_samples_markers_source[["hart"]][[marker]][["high"]], # vs. high Hart
                           fc_x_group_all_samples_markers_source[["icgc"]][[marker]][["low"]], # low ICGC 
                           fc_x_group_all_samples_markers_source[["hart"]][[marker]][["low"]] # vs. low Hart
                           )
  colnames(df_to_plot) <- c("fc_enh_icgc", "fc_enh_hart", "fc_high_icgc", "fc_high_hart", "fc_low_icgc", "fc_low_hart")
  df_to_plot_long <- df_to_plot %>% pivot_longer(everything(), names_to = "seq_type", values_to = "mut_rate")
  df_to_plot_long$seq_type <- factor(df_to_plot_long$seq_type, levels = c("fc_enh_icgc", "fc_enh_hart", "fc_high_icgc", "fc_high_hart", "fc_low_icgc", "fc_low_hart"))
  
  p <- df_to_plot_long %>% 
    ggplot(., aes(x = seq_type, y = mut_rate, group = seq_type))+
    geom_jitter(width = 0.3, alpha = 0.4)+
    geom_boxplot(alpha=0.2, width = 0.4, outliers = F)+
    stat_compare_means(label = "p.signif", method = "wilcox.test", size = 3, 
                       comparisons = list(c("fc_enh_icgc", "fc_enh_hart"), c("fc_high_icgc", "fc_high_hart"), c("fc_low_icgc", "fc_low_hart")))+
    theme_light()+
    geom_hline(yintercept = 0, col = "grey", linetype = 2)+
    labs(x = "", y = "log2(Fold-change over background)", 
         title = paste0("Mutation frequency in ", marker, " enhancers"), 
         subtitle = "ICGC vs. Hartwig")+
    stat_summary(
      fun = mean, 
      geom = "text", 
      aes(label = sprintf("%.2f", ..y..)), 
      vjust = -20,
      color = "blue"
    )
  print(p)
  
  # Same, but log
  df_to_plot_log <- cbind.fill(log2(fc_all_samples_markers_source[["icgc"]][[marker]]+0.1), 
                           log2(fc_all_samples_markers_source[["hart"]][[marker]]+0.1),
                           log2(fc_x_group_all_samples_markers_source[["icgc"]][[marker]][["high"]]+0.1), # high ICGC 
                           log2(fc_x_group_all_samples_markers_source[["hart"]][[marker]][["high"]]+0.1), # vs. high Hart
                           log2(fc_x_group_all_samples_markers_source[["icgc"]][[marker]][["low"]]+0.1), # low ICGC 
                           log2(fc_x_group_all_samples_markers_source[["hart"]][[marker]][["low"]]+0.1) # vs. low Hart
  )
  colnames(df_to_plot_log) <- c("fc_enh_icgc", "fc_enh_hart", "fc_high_icgc", "fc_high_hart", "fc_low_icgc", "fc_low_hart")
  df_to_plot_log_long <- df_to_plot_log %>% pivot_longer(everything(), names_to = "seq_type", values_to = "mut_rate")
  df_to_plot_log_long$seq_type <- factor(df_to_plot_log_long$seq_type, levels = c("fc_enh_icgc", "fc_enh_hart", "fc_high_icgc", "fc_high_hart", "fc_low_icgc", "fc_low_hart"))
  
  p <- df_to_plot_log_long %>% 
    ggplot(., aes(x = seq_type, y = mut_rate, group = seq_type))+
    geom_jitter(width = 0.3, alpha = 0.4)+
    geom_boxplot(alpha=0.2, width = 0.4, outliers = F)+
    stat_compare_means(label = "p.signif", method = "wilcox.test", size = 3, 
                       comparisons = list(c("fc_enh_icgc", "fc_enh_hart"), c("fc_high_icgc", "fc_high_hart"), c("fc_low_icgc", "fc_low_hart")))+
    theme_light()+
    geom_hline(yintercept = 0, col = "grey", linetype = 2)+
    labs(x = "", y = "log2(Fold-change over background)", 
         title = paste0("Mutation frequency in ", marker, " enhancers"), 
         subtitle = "ICGC vs. Hartwig")+
    stat_summary(
      fun = mean, 
      geom = "text", 
      aes(label = sprintf("%.2f", ..y..)), 
      vjust = -20,
      color = "blue"
    )
  print(p)
  
  
}


##

dev.off()

##



