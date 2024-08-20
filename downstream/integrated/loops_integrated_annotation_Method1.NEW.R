
### loops_integrated_annotation_Method1 ###
#'
#' @description
#' This script is used to annotate loops to GRHL2 enhancers and DEGs promoters 
#' Loops were annotated according to 'Method 1', namely considering a window of +/- 1 bin
#' to define if a loop bin overlaps one enhancer or one DEG promoter

library(fs)  
library(tidyverse)
library(GenomicRanges)
library(viridis)
library(gridExtra)

SEED <- 4321
set.seed(SEED)
location <- "local" # 'local' or 'hpc'

### Hi-ChIP Loops ###
kb <- 2

source("/Users/ieo6983/Desktop/fragile_enhancer_clinical/utils/loops_functions.R")
path_hichip <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/HiChip/filtered_loops/")
path_enhancers <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/downstream/")
path_enhancers_ctip <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/downstream/CtIP_enh.hq_signal.clustered.tsv")
path_enhancers_grhl <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/downstream/GRHL_enh.hq_signal.clustered.tsv")

path_tss <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/others/TSSs_elisa/TSSs_from_USCS_hg19_EMSEMBL.tsv")
path_degs <- fs::path("/Users/ieo6983/Desktop/expression/DEGs/Df_DEGs.df_LFC_sig.padj_0.05.log2FC_1.Up_and_Down.tsv")

path_results <- fs::path(paste0("/Users/ieo6983/Desktop/fragile_enhancer_clinical/results/integrated/NEW/", kb, "kb"))

if(location == "hpc"){
  source("/hpcnfs/scratch/PGP/Ciacci_et_al/fragile_enhancer_clinical/utils/loops_functions.R")
  path_hichip <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/HiChip/filtered_loops/")
  path_enhancers <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/downstream/")
  path_enhancers_ctip <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/downstream/CtIP_enh.hq_signal.clustered.tsv")
  path_enhancers_grhl <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/downstream/GRHL_enh.hq_signal.clustered.tsv")
 
  path_tss <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/others/TSSs_elisa/TSSs_from_USCS_hg19_EMSEMBL.tsv")
  path_degs <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/results/expression/RNA_seq/DEGs/Df_DEGs.df_LFC_sig.padj_0.05.log2FC_1.Up_and_Down.tsv")
  
  path_results <- fs::path(paste0("/hpcnfs/scratch/PGP/Ciacci_et_al/results/integrated/NEW/", kb, "kb"))
}


##


# Read loops
scr <- read_tsv(fs::path(path_hichip, 'SCR', sprintf('SCR_loops_%s_kb.tsv', kb)))
kd <- read_tsv(fs::path(path_hichip, 'KD', sprintf('KD_loops_%s_kb.tsv', kb)))
dim(scr); dim(kd)


# SCR
scr_specific <- anti_join(scr, kd, by = c('seqnames1', 'start1', 'end1', 'seqnames2', 'start2', 'end2'))

# KD
kd_specific <- anti_join(kd, scr, by = c('seqnames1', 'start1', 'end1', 'seqnames2', 'start2', 'end2'))

# Shared
common <- inner_join(scr, kd, by = c('seqnames1', 'start1', 'end1', 'seqnames2', 'start2', 'end2'), 
                     suffix = c("_SCR", "_KD"))

# checks
dim(scr_specific)[1] + dim(common)[1] == dim(scr)[1] 
dim(kd_specific)[1] + dim(common)[1] == dim(kd)[1]
sum(duplicated(scr)); sum(duplicated(kd)); sum(duplicated(common));

# Save condition-specific loops
#scr_specific %>% write_tsv(., fs::path(path_results, paste0("/data/", kb, "kb_SCR_specific_loops.tsv")))
#kd_specific %>% write_tsv(., fs::path(path_results, paste0("/data/", kb, "kb_KD_specific_loops.tsv")))


##


## Loops involving GRHL2-enhancers  

# Read GRHL2-bound enhancers
enh_grhl2 <- read_tsv(path_enhancers_grhl)

# change enhancers name  
enh_grhl2$chrom <- paste0("chr", enh_grhl2$chrom)
enh_grhl2$name <- str_c(enh_grhl2$chrom, enh_grhl2$summit, sep=":")

# SCR
scr_enh_overlaps <- find_loops_overlapping_enhancers(enhancers_file = enh_grhl2, loops_file = scr, 
                                                bin_len_kb = kb, ext_nbins = 1)

# KD 
kd_enh_overlaps <- find_loops_overlapping_enhancers(enhancers_file = enh_grhl2, loops_file = kd, 
                                                     bin_len_kb = kb, ext_nbins = 1)

# Loops involving at least 1 GRHL2-bound enhancer (either bin1 or bin2)
scr_enh <- scr_enh_overlaps %>% dplyr::filter(., !is.na(name1) | !is.na(name2))
kd_enh <- kd_enh_overlaps %>% dplyr::filter(., !is.na(name1) | !is.na(name2))

scr_enh_specific <- inner_join(scr_enh, scr_specific, by = c('seqnames1', 'start1', 'end1', 'seqnames2', 'start2', 'end2'))
kd_enh_specific <- inner_join(kd_enh, kd_specific, by = c('seqnames1', 'start1', 'end1', 'seqnames2', 'start2', 'end2'))

common_enh <- inner_join(scr_enh, kd_enh, by = c('seqnames1', 'start1', 'end1', 'seqnames2', 'start2', 'end2'), 
                     suffix = c("_SCR", "_KD"))

dim(scr_enh)[1] == dim(scr_enh_specific)[1] + dim(common_enh)[1]
dim(kd_enh)[1] == dim(kd_enh_specific)[1] + dim(common_enh)[1]


##


## Some stats

# Number of loops per group
df1 <- data.frame("var" = c("scr_tot", "scr_tot_specific", "kd_tot", "kd_tot_specific", "common"))
df1$var <- factor(df1$var, levels = c("scr_tot", "scr_tot_specific", "kd_tot", "kd_tot_specific", "common"))
df1$value <- c(dim(scr)[1], dim(scr_specific)[1], dim(kd)[1], dim(kd_specific)[1], dim(common)[1])
df1$group <- c(rep("scr", 2), rep("kd", 2), "common")
df1$group <- factor(df1$group, levels = c("scr", "kd", "common"))

p1 <- df1 %>% ggplot(., aes(x = group, y = value, fill = var, order = group, width = 0.6))+
  geom_bar(stat = "identity", position = "dodge")+
  theme_light()+
  labs(title = "Number of loops per group", 
       y = "", x = "")+
  scale_y_continuous(n.breaks = 10)+
  scale_fill_viridis(discrete = T)
print(p1)


# Number of GRHL2-bound enhancers loops x group
df2 <- data.frame("var" = c("scr_enh", "scr_enh_specific", "kd_enh", "kd_enh_specific", "common_enh"))
df2$var <- factor(df2$var, levels = c("scr_enh", "scr_enh_specific", "kd_enh", "kd_enh_specific", "common_enh"))
df2$value <- c(dim(scr_enh)[1], dim(scr_enh_specific)[1], dim(kd_enh)[1], dim(kd_enh_specific)[1], dim(common_enh)[1])
df2$group <- c(rep("scr_enh", 2), rep("kd_enh", 2), "common_enh")
df2$group <- factor(df2$group, levels = c("scr_enh", "kd_enh", "common_enh"))

df2$perc <- round(df2$value / df1$value * 100)

p2 <- df2 %>% ggplot(., aes(x = group, y = value, fill = var, order = group, 
                            width = 0.6, label = paste0(perc, "%")))+
  geom_bar(stat = "identity", position = "dodge")+
  geom_text(position = position_dodge(width = 0.65), vjust = -0.2)+
  theme_light()+
  labs(title = "Number of GRHL2-enhancers loops per group", 
       y = "", x = "")+
  scale_y_continuous(n.breaks = 10)+
  scale_fill_viridis(discrete = T)
print(p2)


# Alltogether
df3 <- rbind(df1, df2[,-4])
df3$group2 <- rep(c("scr", "scr_specific", "kd", "kd_specific", "common"), 2)
df3$group2 <- factor(df3$group2, levels = c("scr", "scr_specific", "kd", "kd_specific", "common"))

p3 <- df3 %>% ggplot(., aes(x = group2, y = value, fill = var, order = group, width = 0.6))+
  geom_bar(stat = "identity", position = "dodge")+
  theme_light()+
  labs(title = "Number of GRHL2-enhancers loops per group", 
       y = "", x = "")+
  scale_y_continuous(n.breaks = 10)+
  scale_fill_viridis(discrete = T)
  #facet_wrap(~group)
print(p3)

#ggsave(plot = p1, filename = fs::path(path_results, sprintf("/plots/%skb_number_of_loops",kb), ext = "png"), device = "png", width = 7, height = 5)
#ggsave(plot = p2, filename = fs::path(path_results, sprintf("/plots/%skb_number_of_GRHL2_enhancers_loops",kb), ext = "png"), device = "png", width = 7, height = 5)
#ggsave(plot = p3, filename = fs::path(path_results, sprintf("/plots/%skb_number_of_GRHL2_enhancers_loops_over_tot",kb), ext = "png"), device = "png", width = 7, height = 5)

df2$group <- factor(df2$group, levels = rev(levels(df2$group)))
df2$var <- factor(df2$var, levels = rev(levels(df2$var)))
p <- df2 %>% ggplot(., aes(x = group, y = value, fill = var, order = group, width = 0.6))+
  geom_bar(stat = "identity", position = "dodge")+
  theme_light()+
  labs(title = "Number of GRHL2-enhancers loops per group", 
       y = "", x = "")+
  scale_fill_viridis(discrete = T, direction = -1)+
  coord_flip()
p
#ggsave(plot = p, filename = fs::path(path_results, sprintf("/plots/%skb_number_of_GRHL2_enhancers_loops.flip",kb), ext = "png"), device = "png", width = 7, height = 5)


##


## DEGs Promoters overlap 

# Read TSSs
TSSs <- read_tsv(path_tss)
DEGs <- read_tsv(path_degs)

# All DEGs are present in TSS file
sum(!toupper(DEGs$gene_id) %in% toupper(TSSs$gene_id))

# Add TSSs to DEGs
DEGs_tss <- left_join(DEGs, TSSs[,c(1:2,4)], by = "gene_id")

# Multiple TSSs for each gene
n_tss <- DEGs_tss %>% group_by(gene_id) %>%
  summarise(n_tss = length(unique(tss)))
table(n_tss$n_tss)
# Keep all. When assigninng to loops, tss information is lost and only gene_name remains. Remove duplicates later.

# No ambiguous gene names
n_ids <- DEGs_tss %>% group_by(gene_id) %>%
  summarise(n_ids = length(unique(gene_name)))
table(n_ids$n_ids)

# Find promoter overlaps
scr_enh_degs <- find_additional_DEGS_overlap(DEGs_tss_file = DEGs_tss, loops_enh_file = scr_enh, bin_len_kb = kb)
kd_enh_degs <- find_additional_DEGS_overlap(DEGs_tss_file = DEGs_tss, loops_enh_file = kd_enh, bin_len_kb = kb)

sum(duplicated(scr_enh_degs))
sum(duplicated(kd_enh_degs))

# How many loops have ambiguous annotations in at least 1 bin? 
perc_amb1 <- dim(scr_enh_degs[!is.na(scr_enh_degs$name1) & !is.na(scr_enh_degs$gene_name1) | !is.na(scr_enh_degs$name2) & !is.na(scr_enh_degs$gene_name2), ])[1]/
  dim(scr_enh_degs)[1] *100
perc_amb2 <- dim(kd_enh_degs[!is.na(kd_enh_degs$name1) & !is.na(kd_enh_degs$gene_name1) | !is.na(kd_enh_degs$name2) & !is.na(kd_enh_degs$gene_name2), ])[1]/
  dim(kd_enh_degs)[1] *100

print(sprintf("~%d%% of GRHL2 enhancer-associated loops have at least 1 bin with ambiguous annotation in SCR", round(perc_amb1)))
print(sprintf("~%d%% of GRHL2 enhancer-associated loops have at least 1 bin with ambiguous annotation in KD", round(perc_amb2)))


##


# Save tables
#scr_enh_degs %>% write_tsv(., fs::path(path_output, sprintf("/data/%skb_SCR.anno_loops.GRHL2_enh_DEGs_prom",kb), ext = "tsv"))
#kd_enh_degs %>% write_tsv(., fs::path(path_output, sprintf("/data/%skb_KD.anno_loops.GRHL2_enh_DEGs_prom",kb), ext = "tsv"))


##


# Some stats 

scr_enh_degs_spec <- inner_join(scr_enh_degs, scr_specific, by = c('seqnames1', 'start1', 'end1', 'seqnames2', 'start2', 'end2'))
kd_enh_degs_spec <- inner_join(kd_enh_degs, kd_specific, by = c('seqnames1', 'start1', 'end1', 'seqnames2', 'start2', 'end2'))
enh_degs_common <- inner_join(scr_enh_degs, kd_enh_degs, by = c('seqnames1', 'start1', 'end1', 'seqnames2', 'start2', 'end2'))

dim(scr_enh_degs)[1] == dim(scr_enh_degs_spec)[1]+dim(enh_degs_common)[1]
dim(kd_enh_degs)[1] == dim(kd_enh_degs_spec)[1]+dim(enh_degs_common)[1]

# Loops involving one GRHL2-enhancer and 1 DEG
scr_enh_degs_only <- scr_enh_degs %>% filter(!((is.na(name1) & is.na(gene_name1)) | (is.na(name2) & is.na(gene_name2)))) %>% 
  filter(!is.na(gene_name1) | !is.na(gene_name2)) 
kd_enh_degs_only <- kd_enh_degs %>% filter(!((is.na(name1) & is.na(gene_name1)) | (is.na(name2) & is.na(gene_name2)))) %>% 
  filter(!is.na(gene_name1) | !is.na(gene_name2)) 

scr_enh_degs_only_spec <- inner_join(scr_enh_degs_only , scr_specific, by = c('seqnames1', 'start1', 'end1', 'seqnames2', 'start2', 'end2'))
kd_enh_degs_only_spec <- inner_join(kd_enh_degs_only , kd_specific, by = c('seqnames1', 'start1', 'end1', 'seqnames2', 'start2', 'end2'))

enh_degs_only_common <- inner_join(scr_enh_degs_only, kd_enh_degs_only, by = c('seqnames1', 'start1', 'end1', 'seqnames2', 'start2', 'end2'))

dim(scr_enh_degs_only)[1] == dim(scr_enh_degs_only_spec)[1]+dim(enh_degs_only_common)[1]
dim(kd_enh_degs_only)[1] == dim(kd_enh_degs_only_spec)[1]+dim(enh_degs_only_common)[1]

# Save tables
#scr_enh_degs_only %>% write_tsv(., fs::path(path_results, paste0("/data/", kb, "kb_SCR.anno_loops.GRHL2_enh_DEGs_prom_ONLY.tsv")))
#kd_enh_degs_only %>% write_tsv(., fs::path(path_results, paste0("/data/", kb, "kb_KD.anno_loops.GRHL2_enh_DEGs_prom_ONLY.tsv")))

# Number of enh-prom loops x group
dfs <- list("scr" = list("tot" = scr_enh_degs, "spec" = scr_enh_degs_spec, 
                         "only" = scr_enh_degs_only, "only_spec" = scr_enh_degs_only_spec), 
            "kd" = list("tot" = kd_enh_degs, "spec" = kd_enh_degs_spec, 
                        "only" = kd_enh_degs_only, "only_spec" = kd_enh_degs_only_spec))
df4 <- list("scr" = df4_scr <- as.data.frame(matrix(nrow = 1, ncol = 1)), 
            "kd" = df4_kd <- as.data.frame(matrix(nrow = 1, ncol = 1)))

for(cond in c("scr", "kd")){
  dfs_cond <- dfs[[cond]]
  
  #df4[[cond]]$group2 <- factor(c(paste0(cond,"_tot"), paste0(cond,"_spec")), levels = c(paste0(cond,"_tot"), paste0(cond,"_spec")))

  # Number of annotated loops
  df4[[cond]]$anno_loops <- dim(dfs[[cond]][["tot"]])[1]
  df4[[cond]]$anno_loops_cond_spec <- dim(dfs[[cond]][["spec"]])[1]
  
  # Number of ENH-DEG loops only
  df4[[cond]]$enh_deg_loops <- dim(dfs[[cond]][["only"]])[1]
  df4[[cond]]$enh_deg_loops_spec <- dim(dfs[[cond]][["only_spec"]])[1]
  
  # Number of GRHL2 enhancers - both bins
  df4[[cond]]$enh_n_uniq <- length(unique(c(na.omit(dfs[[cond]][["only"]]$name1), na.omit(dfs[[cond]][["only"]]$name2))))
  # Number of DEGs
  df4[[cond]]$degs_n_uniq <- length(unique(c(na.omit(dfs[[cond]][["only"]]$gene_name1), na.omit(dfs[[cond]][["only"]]$gene_name2))))
  
  # Number of Up / Down DEGs
  up1 <- dfs[[cond]][["only"]] %>% filter(DE1 == "Up")
  up2 <- dfs[[cond]][["only"]] %>% filter(DE2 == "Up")
  df4[[cond]]$degs_up <- length(unique(na.omit(c(up1$gene_name1, up2$gene_name2))))
  down1 <- dfs[[cond]][["only"]] %>% filter(DE1 == "Down")
  down2 <- dfs[[cond]][["only"]] %>% filter(DE2 == "Down")
  df4[[cond]]$degs_down <- length(unique(na.omit(c(down1$gene_name1, down2$gene_name2))))

}

# Merge scr and kd
df4[["scr"]]$V1 <- "scr"
df4[["kd"]]$V1 <- "kd"
df4 <- rbind(df4[["scr"]], df4[["kd"]]) %>% rename(V1 = "group1")

df4
    

##


# Loop types - ALL annotated
dfs <- list("scr" = list("tot" = scr_enh_degs, "only" = scr_enh_degs_only), 
            "kd" = list("tot" = kd_enh_degs, "only" = kd_enh_degs_only))
df5 <- list("scr" = df5_scr <- as.data.frame(matrix(nrow = 1, ncol = 1)), 
            "kd" = df5_kd <- as.data.frame(matrix(nrow = 1, ncol = 1)))

for(cond in c("scr", "kd")){
  df_cond <- dfs[[cond]]
  
  df5[[cond]]$tot <- dim(df_cond[["tot"]])[1]
  
  # Loops connecting 1 enhancer with 1 DEG promoter
  df5[[cond]]$enh_deg <- dim(df_cond[["only"]])[1]
  
  # Loops connecting 1 enhancer with 1 enhancer
  enh_enh <- df_cond[["tot"]] %>% filter((!is.na(name1) & !is.na(name2))) 
                    #%>% dplyr::filter((is.na(gene_name1) & is.na(gene_name2)))
  df5[[cond]]$enh_enh <- dim(enh_enh)[1]
  
  # Loops connecting 1 DEG promoter with 1 DEG promoter
  prom_prom <- df_cond[["tot"]] %>% filter((!is.na(gene_name1) & !is.na(gene_name2)))
                #%>% dplyr::filter((is.na(name1) & is.na(name2)))
  df5[[cond]]$prom_prom <- dim(prom_prom)[1]
  
  # Loops connecting 1 enhancer with "any" (not-annotated region)
  enh_any <- df_cond[["tot"]] %>% filter((!is.na(name1) & is.na(name2) & is.na(gene_name2)) | (!is.na(name2) & is.na(name1) & is.na(gene_name1)))
  df5[[cond]]$enh_any <- dim(enh_any)[1]
  
  # Loops connecting 1 DEG promoter with "any" (noot-annotated region)
  prom_any <- df_cond[["tot"]] %>% filter((!is.na(gene_name1) & is.na(gene_name2) & is.na(name2)) | (!is.na(gene_name2) & is.na(gene_name1) & is.na(name1)))
  df5[[cond]]$prom_any <- dim(prom_any)[1]
  
  # Ambiguous loops: at least one bin overlaps with both enh and promoter
  amb <- df_cond[["tot"]] %>% filter(!is.na(name1) & !is.na(gene_name1) | !is.na(name2) & !is.na(gene_name2))
  df5[[cond]]$amb <- dim(amb)[1]
}

# Merge scr and kd
df5[["scr"]]$V1 <- "scr"
df5[["kd"]]$V1 <- "kd"
df5 <- rbind(df5[["scr"]], df5[["kd"]]) %>% rename(V1 = "group1")

df5


##


# Loops enh-degs only
dfs <- list("scr" = scr_enh_degs_only, "kd" = kd_enh_degs_only)
df6 <- list("scr" = df6_scr <- as.data.frame(matrix(nrow = 1, ncol = 1)), 
            "kd" = df6_kd <- as.data.frame(matrix(nrow = 1, ncol = 1)))

for(cond in c("scr", "kd")){
  df_cond <- dfs[[cond]]
  
  df6[[cond]]$tot <- dim(df_cond)[1]
    
  # Loops connecting 1 enhancer with 1 DEG promoter
  enh_deg <- df_cond %>% filter((!is.na(name1) & is.na(gene_name1)) & (is.na(name2) & !is.na(gene_name2)) | (!is.na(name2) & is.na(gene_name2)) & (is.na(name1) & !is.na(gene_name1)))  
  df6[[cond]]$enh_deg <- dim(enh_deg)[1]
  
  # Loops connecting 1 enhancer with 1 enhancer
  enh_enh <- df_cond %>% filter((!is.na(name1) & !is.na(name2)) & (is.na(gene_name1) & is.na(gene_name2)))
  df6[[cond]]$enh_enh <- dim(enh_enh)[1]
  
  # Loops connecting 1 DEG promoter with 1 DEG promoter
  prom_prom <- df_cond %>% filter((is.na(name1) & is.na(name2)) & (!is.na(gene_name1) & !is.na(gene_name2)))
  df6[[cond]]$prom_prom <- dim(prom_prom)[1]
  
  # Loops connecting 1 enhancer with "any" (not-annotated region)
  enh_any <- df_cond %>% filter((!is.na(name1) & is.na(gene_name1) & is.na(name2) & is.na(gene_name2))
                                | (!is.na(name2) & is.na(gene_name2) & is.na(name1) & is.na(gene_name1)))
  df6[[cond]]$enh_any <- dim(enh_any)[1]
    
  # Loops connecting 1 DEG promoter with "any" (noot-annotated region)
  prom_any <- df_cond %>% filter((!is.na(gene_name1) & is.na(name1) & is.na(gene_name2) & is.na(name2)) 
                                 | !is.na(gene_name2) & is.na(name2) & is.na(gene_name1) & is.na(name1))
  df6[[cond]]$prom_any <- dim(prom_any)[1]
  
  # Ambiguous loops: at least one bin overlaps with both enh and promoter
  amb <- df_cond %>% filter(!is.na(df_cond$name1) & !is.na(df_cond$gene_name1) | !is.na(df_cond$name2) & !is.na(df_cond$gene_name2))
  df6[[cond]]$amb <- dim(amb)[1]
}

# Merge scr and kd
df6[["scr"]]$V1 <- "scr"
df6[["kd"]]$V1 <- "kd"
df6 <- rbind(df6[["scr"]], df6[["kd"]]) %>% rename(V1 = "group1")

df6


##


# Extract DEGs

one <- scr_enh_degs_only[, c("gene_name1", "DE1")][!duplicated(scr_enh_degs_only[, c("gene_name1", "DE1")]), ]
colnames(one) <- c("gene_name", "DE")
two <- scr_enh_degs_only[, c("gene_name2", "DE2")][!duplicated(scr_enh_degs_only[, c("gene_name2", "DE2")]), ]
colnames(two) <- c("gene_name", "DE")
full <- rbind(one, two)
full <- full[!duplicated(full), ]
full <- full %>% arrange(., DE)
#write_tsv(full, fs::path(path_results, paste0("DEGs_SCR_", kb, "kb_M1.tsv")))

one <- kd_enh_degs_only[, c("gene_name1", "DE1")][!duplicated(kd_enh_degs_only[, c("gene_name1", "DE1")]), ]
two <- kd_enh_degs_only[, c("gene_name2", "DE2")][!duplicated(kd_enh_degs_only[, c("gene_name2", "DE2")]), ]
colnames(one) <- c("gene_name", "DE"); colnames(two) <- c("gene_name", "DE")
full <- rbind(one, two)
full <- full[!duplicated(full), ]
full <- full %>% arrange(., DE)
#write_tsv(full, fs::path(path_results, paste0("DEGs_KD_", kb, "kb_M1.tsv")))



