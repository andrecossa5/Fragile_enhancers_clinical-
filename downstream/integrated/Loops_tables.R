
#' Annotate loops to enhancers and DEGs promoters
#'
#' @description
#' This script is used to annotate loops to GRHL2 enhancers and DEGs promoters
#' And to generate a unique table containing both SCR and KD loops 
#' 


library(fs)  # File manipulations
library(tidyverse)

SEED <- 4321
set.seed(SEED)
#source("./Desktop/enhancers_project/Analyses/loops/loops_functions.R")


### Loops Table - Method 1 ###
kb <- 4

path_main <- "/Users/ieo6983/Desktop/fragile_enhancer_clinical"
path_data <- fs::path(path_main, "data") 
path_results <- fs::path(path_main, "results")
path_hichip <- fs::path(path_data, "functional_genomics/HiChip/filtered_loops/")
path_enhancers <- fs::path(path_data, "functional_genomics/others")
path_tss <- file.path(path_data, "functional_genomics/others/TSSs_elisa/TSSs_from_USCS_hg19_EMSEMBL.tsv")
path_degs <- "/Users/ieo6983/Desktop/expression/DEGs/Df_DEGs.df_LFC_sig.padj_0.05.log2FC_1.Up_and_Down.tsv"

path_main_input <- sprintf("/Users/ieo6983/Desktop/enhancers_project/Analyses/loops/results/%skb/data", kb)
path_output <- fs::path(sprintf("/Users/ieo6983/Desktop/enhancers_project/Analyses/loops/results/%skb", kb))


##


#' Find and filter loops overlapping with enhancers - Method 1, modified
#'
#' This function annotates loops' bins to a set of input enhancers
#' 
#' @param enhancer_file File with enhancers coordinates: chrom, start, end, summit, name, (cluster). I.e. GRHL2-enhancers
#' @param loops_file File with loops coordinates: 'seqnames1', 'start1', 'end1', 'seqnames2', 'start2', 'end2', others (p-value, FDR, ..)
#' @param bin_len_kb Length of loop bins resulting from loop calling in kb. I.e. 4, for 4 kb.
#' @param ext_bins Number of bins left/right by which loop bin is extended. Overlaps are searched within the extended bin.
#' @examples
#' scr_enh_overlaps <- find_loops_overlapping_enhancers(enhancers_file = enh_grhl2, loops_file = scr, bin_len_kb = 4, ext_nbins = 1)

find_loops_overlapping_enhancers <- function(enhancers_file, 
                                             loops_file, 
                                             bin_len_kb = 4,
                                             ext_nbins = 1){
  
  # Conver input files to GRanges objects
  enhancers_gr <- makeGRangesFromDataFrame(enhancers_file, keep.extra.columns = T)
  loops_bin1 <- loops_file[, endsWith(colnames(loops_file), "1")]
  loops_bin2 <- loops_file[, endsWith(colnames(loops_file), "2")]
  cols <- colnames(loops_bin1)
  loops_bin1_gr <- makeGRangesFromDataFrame(loops_bin1, seqnames.field = cols[1], 
                                            start.field = cols[2], end.field = cols[3])
  cols <- colnames(loops_bin2)
  loops_bin2_gr <- makeGRangesFromDataFrame(loops_bin2, seqnames.field = cols[1], 
                                            start.field = cols[2], end.field = cols[3])
  
  # Check if all regions have equal width, equal to kb
  if(all(width(loops_bin1_gr) == width(loops_bin1_gr)[1])){
    message(sprintf("All loop bins 1 have equal length; equal to: %d", width(loops_bin1_gr)[1]))
  } else {
    message("WARNING: not all loop bins 1 have equal length")
  }
  
  if(all(width(loops_bin2_gr) == width(loops_bin2_gr)[1])){
    message(sprintf("All loop bins 2 have equal length; equal to: %d", width(loops_bin2_gr)[1]))
  } else {
    message("WARNING: not all loop bins 1 have equal length")
  }
  
  # Extend loop bins of +/- n bins
  start(loops_bin1_gr) <- start(loops_bin1_gr) - (bin_len_kb*1000*ext_nbins)
  end(loops_bin1_gr) <- end(loops_bin1_gr) + (bin_len_kb*1000*ext_nbins)
  start(loops_bin2_gr) <- start(loops_bin2_gr) - (bin_len_kb*1000*ext_nbins)
  end(loops_bin2_gr) <- end(loops_bin2_gr) + (bin_len_kb*1000*ext_nbins)
  
  # Find overlaps 
  #bin1
  ovrlp_bin1 <- findOverlaps(loops_bin1_gr, enhancers_gr, select = "first")
  enh_1 <- enhancers_file[ovrlp_bin1, c("name", "cluster")]
  colnames(enh_1) <- paste0(colnames(enh_1), 1)
  loops_bin1_ovrlp <- cbind(loops_bin1, enh_1) # GRanges are extended, but df remains as originally
  #bin2
  ovrlp_bin2 <- findOverlaps(loops_bin2_gr, enhancers_gr, select = "first")
  enh_2 <- enhancers_file[ovrlp_bin2, c("name", "cluster")]
  colnames(enh_2) <- paste0(colnames(enh_2), 2)
  loops_bin2_ovrlp <- cbind(loops_bin2, enh_2) # GRanges are extended, but df remains as originally
  #joined
  loops_ovrlp <- cbind(loops_bin1_ovrlp, loops_bin2_ovrlp) %>%
    relocate(., c(name1,cluster1), .after = end2)
  
  # Added to keep counts and FDR information for each loop
  loops_ovrlp <- cbind(loops_ovrlp, loops_file[,c("counts","qvalue")])
  
  return(loops_ovrlp)
}

#' Filter loops overlapping with DEGs, modified
#' 
#' This function ...
#' 
#' @param name description
#' @examples
#' # example code
# DEGs_tss_file:   File with DEGs TSSs coordinates: 
# loops_enh_file:   File with coordinates of loops, with info on overlapping enahncers: 'seqnames1', 'start1', 'end1', 'seqnames2', 'start2', 'end2', name1, cluster1, ...
# bin_len_kb:       Length of loop bins resulting from loop calling in kb. I.e. 4, for 4 kb.
# ext_nbins:        Number of bins left/right by which loop bin is extended. Overlaps are searched within the extended bin.

find_additional_DEGS_overlap <- function(DEGs_tss_file, 
                                         loops_enh_file, 
                                         bin_len_kb = 4,
                                         ext_nbins = 1){
  
  tss_file <- DEGs_tss_file
  # Create GRanges obj for TSSs
  tss_file$start <- tss_file$tss
  tss_file$end <- tss_file$tss
  tss_file <- tss_file %>% dplyr::select(gene_name, chrom, tss, start, end, DE, log2FoldChange, padj, pos_tss)
  tss_gr <- makeGRangesFromDataFrame(tss_file, keep.extra.columns = T, seqnames.field = "chrom")
  
  # Create GRanges obj for loops_enh
  loops_bin1 <- loops_enh_file[, endsWith(colnames(loops_enh_file), "1")]
  loops_bin2 <- loops_enh_file[, endsWith(colnames(loops_enh_file), "2")]
  cols <- colnames(loops_bin1)
  loops_bin1_gr <- makeGRangesFromDataFrame(loops_bin1, seqnames.field = cols[1], 
                                            start.field = cols[2], end.field = cols[3])
  cols <- colnames(loops_bin2)
  loops_bin2_gr <- makeGRangesFromDataFrame(loops_bin2, seqnames.field = cols[1], 
                                            start.field = cols[2], end.field = cols[3])
  
  
  # Extend loop bins of +/- n bins
  start(loops_bin1_gr) <- start(loops_bin1_gr) - (bin_len_kb*1000*ext_nbins)
  end(loops_bin1_gr) <- end(loops_bin1_gr) + (bin_len_kb*1000*ext_nbins)
  start(loops_bin2_gr) <- start(loops_bin2_gr) - (bin_len_kb*1000*ext_nbins)
  end(loops_bin2_gr) <- end(loops_bin2_gr) + (bin_len_kb*1000*ext_nbins)
  
  # Find overlaps 
  #bin1
  ovrlp_bin1 <- findOverlaps(loops_bin1_gr, tss_gr, select = "first")
  tss_1 <- tss_file[ovrlp_bin1, c("gene_name", "DE", "log2FoldChange", "padj", "pos_tss")]
  colnames(tss_1) <- paste0(colnames(tss_1), 1)
  loops_bin1_ovrlp <- cbind(loops_bin1, tss_1) # GRanges are extended, but df remains as originally
  #bin2
  ovrlp_bin2 <- findOverlaps(loops_bin2_gr, tss_gr, select = "first")
  tss_2 <- tss_file[ovrlp_bin2, c("gene_name", "DE", "log2FoldChange", "padj", "pos_tss")]
  colnames(tss_2) <- paste0(colnames(tss_2), 2)
  loops_bin2_ovrlp <- cbind(loops_bin2, tss_2) # GRanges are extended, but df remains as originally
  
  #joined
  loops_ovrlp <- cbind(loops_bin1_ovrlp, loops_bin2_ovrlp) %>% 
    relocate(., c(name1,cluster1,gene_name1, DE1, log2FoldChange1, padj1, pos_tss1), .after = end2)
  
  # Added to keep counts and FDR information for each loop
  loops_ovrlp <- cbind(loops_ovrlp, loops_enh_file[,c("counts","qvalue")])

  # Remove eventual duplicates
  message(sprintf("Removing %d duplicated loops", sum(duplicated(loops_ovrlp))))
  loops_ovrlp <- loops_ovrlp[!duplicated(loops_ovrlp), ]
  
  return(loops_ovrlp)
}


#' Unbundle ambigously annotated bins 
#' 
#' This function ...
#' @param name description
#' @examples
#' # example code
  
unbundle_ambiguous_bins <- function(enh_degs_loops_file){
  
  # Count how many overlapping regions x bin
  n1 <- rowSums(!is.na(enh_degs_loops_file[, c("name1", "gene_name1")]))
  cond_amb1 <- n1 > 1
  n2 <- rowSums(!is.na(enh_degs_loops_file[, c("name2", "gene_name2")]))
  cond_amb2 <- n2 > 1
    
  cond_amb <- (n1 > 1 | n2 > 1)
  
  df_not_amb <- enh_degs_loops_file[!cond_amb, ]
  df_amb <- enh_degs_loops_file[cond_amb, ]
  message(sprintf("Check if size of df_amb and df_not_amb correspond to size of whole: %s", 
                  dim(df_not_amb)[1] + dim(df_amb)[1] == dim(enh_degs_loops_file)[1]))
  
  # Store whether bin1 overlaps with 1 or 2 regions (enhancer & promoter)
  df_amb$n_bin1 <- rowSums(!is.na(df_amb %>% dplyr::select(name1, gene_name1)))
  # Store whether bin2 overlaps with 1 or 2 regions (enhancer & promoter)
  df_amb$n_bin2 <- rowSums(!is.na(df_amb %>% dplyr::select(name2, gene_name2)))
  
  # Case 2-0/1: bin1 overlaps with both prom and enh; Case 0/1-2: bin2 overlaps with both prom and enh; Case 2-2: bin1 & 2 overlap with both prom and enh
  df_amb_new <- df_amb
  
  for(i in 1:dim(df_amb)[1]){
    r <- df_amb[i, ]
    
    # If case 1, duplicate row separating prom and enh in bin2
    if((r$n_bin1 == 1 & r$n_bin2 == 2) | (r$n_bin1 == 0 & r$n_bin2 == 2)){
      df_amb_new[i, ]$gene_name2 <- NA
      df_amb_new[i, c("DE2", "log2FoldChange2", "padj2", "pos_tss2")] <- NA
      
      r$name2 <- NA
      r$cluster2 <- NA
      df_amb_new <- rbind(df_amb_new, r)
    }
    
    if((r$n_bin1 == 2 & r$n_bin2 == 1) | (r$n_bin1 == 2 & r$n_bin2 == 0)){
      df_amb_new[i, ]$gene_name1 <- NA
      df_amb_new[i, c("DE1", "log2FoldChange1", "padj1", "pos_tss1")] <- NA
      
      r$name1 <- NA
      r$cluster1 <- NA
      df_amb_new <- rbind(df_amb_new, r)
    }
    
    if(r$n_bin1 == 2 & r$n_bin2 == 2){
      r1 <- r; r2 <- r; r3 <- r
      
      df_amb_new[i, ]$gene_name1 <- NA; df_amb_new[i, ]$gene_name2 <- NA
      df_amb_new[i, c("DE1", "log2FoldChange1", "padj1", "pos_tss1")] <- NA; df_amb_new[i, c("DE2", "log2FoldChange2", "padj2", "pos_tss2")] <- NA
      
      r1$gene_name1 <- NA; r1$name2 <- NA
      r1[, c("DE1", "log2FoldChange1", "padj1", "pos_tss1")] <- NA
      r1$cluster2 <- NA
      df_amb_new <- rbind(df_amb_new, r1)
      
      r2$name1 <- NA; r2$gene_name2 <- NA
      r2$cluster1 <- NA
      r2[, c("DE2", "log2FoldChange2", "padj2", "pos_tss2")] <- NA
      df_amb_new <- rbind(df_amb_new, r2)
      
      r3$name1 <- NA; r3$name2 <- NA
      r3$cluster1 <- NA; r3$cluster2 <- NA
      df_amb_new <- rbind(df_amb_new, r3)
    }
  }
  
  # Check if size of df_amb_new matches df_amb with duplicated rows
  check <- dim(df_amb)[1]*2
  check <- check + (dim(df_amb%>%filter(n_bin1==2&n_bin2==2))[1]*2)
  message(sprintf("Size of new df is compatible with input: %s", check == dim(df_amb_new)[1]))
  
  # Make sure that now no ambiguous bins are present
  df_amb_new$n_bin1 <- rowSums(!is.na(df_amb_new %>% dplyr::select(name1, gene_name1)))
  df_amb_new$n_bin2 <- rowSums(!is.na(df_amb_new %>% dplyr::select(name2, gene_name2)))
  message(sprintf("No more ambiguous bins: %s", sum(c(df_amb_new$n_bin1, df_amb_new$n_bin2) > 1)==0))
  
  # Re-join df_not_amb & df_amb
  df_amb_new$n_bin1 <- NULL; df_amb_new$n_bin2 <- NULL
  df_new <- rbind(df_not_amb, df_amb_new)
  
  return(df_new)
}


##


# Read loops
scr <- read_tsv(fs::path(path_hichip, 'SCR', sprintf('SCR_loops_%s_kb.tsv', kb)))
kd <- read_tsv(fs::path(path_hichip, 'KD', sprintf('KD_loops_%s_kb.tsv', kb)))

# SCR
scr_specific <- anti_join(scr, kd, by = c('seqnames1', 'start1', 'end1', 'seqnames2', 'start2', 'end2'))

# KD
kd_specific <- anti_join(kd, scr, by = c('seqnames1', 'start1', 'end1', 'seqnames2', 'start2', 'end2'))

# Common
common <- inner_join(scr, kd, by = c('seqnames1', 'start1', 'end1', 'seqnames2', 'start2', 'end2'))

dim(common); dim(scr_specific); dim(kd_specific)
dim(common)[1]+dim(scr_specific)[1] == dim(scr)[1]
dim(common)[1]+dim(kd_specific)[1] == dim(kd)[1]


##


# Read GRHL2-bound enhancers
columns_names <- c("chrom", "start", "end", "cluster")
enh_grhl2 <- read_tsv(fs::path(path_enhancers, "Cluster_GRHL_Enh_All.txt"), col_names = columns_names, comment = "#")

# Pre-process: add summit & enhancer name 
enh_grhl2$summit <- enh_grhl2$end
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


##


## DEGs Promoters overlap 

# Read TSSs
TSSs <- read_tsv(path_tss)
DEGs <- read_tsv("./Desktop/expression/DEGs/Df_DEGs.df_LFC_sig.padj_0.05.log2FC_1.Up_and_Down.tsv")

# All DEGs are present in TSS file
sum(!toupper(DEGs$gene_id) %in% toupper(TSSs$gene_id))

# Add TSSs to DEGs
DEGs_tss <- left_join(DEGs, TSSs[,c(1:2,4)], by = "gene_id")
DEGs_tss$pos_tss <- str_c(DEGs_tss$chrom, DEGs_tss$tss, sep=":")

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


##


scr_enh_degs_unb <- unbundle_ambiguous_bins(scr_enh_degs)
kd_enh_degs_unb <- unbundle_ambiguous_bins(kd_enh_degs)

# Check new sizes
dim(scr_enh_degs); dim(scr_enh_degs_unb)
dim(kd_enh_degs); dim(kd_enh_degs_unb)

# Save tables
#scr_enh_degs %>% write_tsv(., fs::path(path_output, sprintf("/data/tables/%skb_SCR.all_anno_loops.ENH_DEGs_any", kb), ext = "tsv"))
#scr_enh_degs_unb %>% write_tsv(., fs::path(path_output, sprintf("/data/tables/%skb_SCR.all_anno_loops.ENH_DEGs_any.unbundled", kb), ext = "tsv"))
#kd_enh_degs %>% write_tsv(., fs::path(path_output, sprintf("/data/tables/%skb_KD.all_anno_loops.ENH_DEGs_any", kb), ext = "tsv"))
#kd_enh_degs_unb %>% write_tsv(., fs::path(path_output, sprintf("/data/tables/%skb_KD.all_anno_loops.ENH_DEGs_any.unbundled", kb), ext = "tsv"))


##


# Unique table 
scr_enh_degs <- read_tsv(fs::path(path_output, sprintf("/data/tables/%skb_SCR.all_anno_loops.ENH_DEGs_any", kb), ext = "tsv"))
kd_enh_degs <- read_tsv(fs::path(path_output, sprintf("/data/tables/%skb_KD.all_anno_loops.ENH_DEGs_any", kb), ext = "tsv"))

scr_enh_degs_spec <- inner_join(scr_enh_degs, scr_specific, by = c('seqnames1', 'start1', 'end1', 'seqnames2', 'start2', 'end2'))
kd_enh_degs_spec <- inner_join(kd_enh_degs, kd_specific, by = c('seqnames1', 'start1', 'end1', 'seqnames2', 'start2', 'end2'))
enh_degs_common <- inner_join(scr_enh_degs, kd_enh_degs, by = c('seqnames1', 'start1', 'end1', 'seqnames2', 'start2', 'end2'))
  
dim(scr_enh_degs_spec); dim(kd_enh_degs_spec); dim(enh_degs_common)
dim(scr_enh_degs_spec)[1] + dim(enh_degs_common)[1] == dim(scr_enh_degs)[1]
dim(kd_enh_degs_spec)[1] + dim(enh_degs_common)[1] == dim(kd_enh_degs)[1]

# Merge the two sets of enhancers into one table only
uni_enh_degs <- full_join(scr_enh_degs, kd_enh_degs, 
                          by = colnames(scr_enh_degs)[1:20], 
                          suffix = c(".scr", ".kd"))

# Do sizes correspond?
dim(scr_enh_degs_spec)[1] + dim(kd_enh_degs_spec)[1] + dim(enh_degs_common)[1] == dim(uni_enh_degs)[1]
#uni_enh_degs %>% write_tsv(., fs::path(path_output, sprintf("/data/tables/%skb_Unified_table.SCR_plus_KD_counts.all_anno_loops.ENH_DEGs_any", kb), ext = "tsv"))


##

