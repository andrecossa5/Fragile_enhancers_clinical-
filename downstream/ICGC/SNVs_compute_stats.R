
#' Stats: 
#' - N. of mutations x enhancers 
#' - N. of donors mutated etc
#' - 




library(tidyverse)
library(GenomicRanges)
library(reshape2)

SEED <- 4321
set.seed(SEED)

source("./utils/functions_genomics.R")

path_SSMs <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/genomics/pre_processed_ICGC/simple_somatic_mutation.open.matching_calls.with_AFs.tsv")
path_enhancers_ctip <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/others/Cluster_CtIP_Enh_All.txt")
path_enhancers_grhl <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/others/Cluster_GRHL_Enh_All.txt")

path_results_data <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/results/ICGC/OLD/enhancers_SSMs_overlaps/data/")
path_results_plots <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/results/ICGC/OLD/enhancers_SSMs_overlaps/plots/")

WIN <- 3000
MARKERS <- c("CtIP", "GRHL") 


##


## Load input files 
SSMs <- read_tsv(path_SSMs)

columns_names <- c("chrom", "start", "end", "cluster")
enh_ctip <- read_tsv(path_enhancers_ctip, col_names = columns_names, comment = "#") 
enh_grhl <- read_tsv(path_enhancers_grhl, col_names = columns_names, comment = "#") 
enh_all <- list("CtIP" = enh_ctip, "GRHL" = enh_grhl)


## Pre-process 

# add summit & name to enhancers 
enh_all <- lapply(enh_all, function(df){
  df$summit <- df$end
  df$name <- str_c(df$chrom, df$summit, sep=":")
  return(df)})

# Filter SSMs (only SNVs)
SNVs <- SSMs %>% rename(., c("chromosome_start" = "start", "chromosome_end" = "end")) %>% 
  mutate(., chromosome = paste0("chr", .$chromosome)) %>%
  dplyr::filter(., end - start == 0)


##


## Find enhancers / SNVs overlaps

input_variants <- SNVs

stat_enhancers_enh <- matrix(nrow = 5000, ncol = length(MARKERS))
colnames(stat_enhancers_enh) <- MARKERS
stat_donors_enh <- matrix(nrow = length(unique(input_variants$icgc_donor_id)), ncol = length(MARKERS))
colnames(stat_donors_enh) <- MARKERS

for(marker in MARKERS){
  print(m)
  marker <- m
  mut <- "SNVs"
  
  ## Pre/process input data 
  marker_enh <- enh_all[[marker]]
  
  # Extend regions of 'win' from summit
  marker_enh$start <- marker_enh$summit - WIN
  marker_enh$end <- marker_enh$summit + WIN
  
  # Generate random input sequences 
  ran_seqs <- generate_random_seqs(SEED, dim(marker_enh)[1], win = WIN)
  
  # Create GRanges objects 
  enh_gr <- makeGRangesFromDataFrame(marker_enh, keep.extra.columns = T)
  ran_seqs_gr = makeGRangesFromDataFrame(ran_seqs, keep.extra.columns = T)
  SNVs_gr <- makeGRangesFromDataFrame(input_variants, keep.extra.columns = T)
  
  
  ##
  
  
  ## Find mutations / enhancers overlaps 
  
  # CtIP/GRHL enhancers 
  hits_obj_enh <- findOverlaps(query=enh_gr, subject=SSMs_gr)
  SSMs_sbj <- data.frame(SSMs_gr[subjectHits(hits_obj_enh)])
  colnames(SSMs_sbj)[1:5] <- paste(colnames(SSMs_sbj)[1:5], "sbj", sep = "_")
  enh_SSMs <- cbind(data.frame(enh_gr[queryHits(hits_obj_enh)]), SSMs_sbj) 
  
  # Random sequences
  hits_obj_ran <- findOverlaps(query=ran_seqs_gr, subject=SSMs_gr)
  SSMs_sbj <- data.frame(SSMs_gr[subjectHits(hits_obj_ran)])
  colnames(SSMs_sbj)[1:5] <- paste(colnames(SSMs_sbj)[1:5], "sbj", sep = "_")
  ran_seqs_SSMs <- cbind(data.frame(ran_seqs_gr[queryHits(hits_obj_ran)], SSMs_sbj))
  
  
  
  ### SAVE stats - Enhancers only
  save_intermediate <- FALSE
  if(save_intermediate == TRUE){
    
    # Save enh_SSMs objects - intermediate object useful for all the ubsequent analyses
    enh_SSMs <- enh_SSMs %>% relocate(., summit, .before = start) 
    enh_SSMs <- enh_SSMs[, -c(8:14)]
    # add Id, useful to compare with RegulomeDB scores
    enh_SSMs$ssm_id <- paste(paste(enh_SSMs$seqnames_sbj, enh_SSMs$start_sbj-1, sep=":"), enh_SSMs$end_sbj,sep="-")
    enh_SSMs <- enh_SSMs %>% relocate(., ssm_id, .before = seqnames_sbj)
    write_tsv(enh_SSMs, paste(OUT_FOLDER_DATA, "Table_enh_SSMs_", marker, ".all_overlaps.tsv",sep=""))
    
    # Save n_enhancers_hit for each donor
    enh_SSMs %>% group_by(icgc_donor_id) %>%
      dplyr::summarise("n_enhancers_hit" = length(unique(name))) %>%
      arrange(., desc(n_enhancers_hit)) %>%
      write_tsv(., paste(OUT_FOLDER_DATA, marker, 
                         "_all_enhancers_SSMs.n_enhancers_hit_x_donor.tsv", sep=""))
    # Save n_donors_hit for each enhancer
    enh_SSMs %>% dplyr::group_by(name) %>%
      dplyr::summarise(., n_donors_hit = length(unique(icgc_donor_id))) %>%
      arrange(., desc(n_donors_hit)) %>%
      write_tsv(., paste(OUT_FOLDER_DATA, marker, "_all_enhancers_SSMs.n_donors_hit.tsv",sep=""))
    
  }
  
  
  
  
###



### Compute stats for plots - all datasets
classes = c(paste(m,"enhancers",sep="_"), "put_enhancers", "ran_seqs")
stat_enh <- compute_stat_enhancers(list(enh_SSMs, put_enh_SSMs, ran_seqs_SSMs), 
                                   list(marker_enh, put_enh, ran_seqs), classes = classes)
stat_don <- compute_stat_donors(list(enh_SSMs, put_enh_SSMs, ran_seqs_SSMs), SSMs, classes = classes)

# Add enhancers stats 
stat_enhancers_enh[, marker] <- c(stat_enh[, classes[1]], 
                                  rep(NA, dim(stat_enhancers_enh)[1] - length(stat_enh[, classes[1]])))
stat_donors_enh[, marker] <- stat_don[, classes[1]]


### Plot stats
avg <- str_flatten(paste(names(apply(stat_enh, 2, mean)), ":", round(apply(stat_enh, 2, mean),2), sep = ""), 
                   collapse = "  ")
perc <- apply(stat_enh, 2, function(x) round(sum(x != 0) / dim(marker_enh)[1] * 100))
perc <- str_flatten(paste(names(perc), ":", perc, " %", sep=""), collapse = "  ")

b1 <- stat_enh %>%
  pivot_longer(everything(), names_to = "cluster", values_to = "value") %>%
  ggplot(., aes(x = cluster, y = value, fill = cluster))+
  geom_boxplot()+
  labs(title = "Number of donors hit", 
       subtitle = paste(
         paste("Average number of donors that each enhancer hits: ", avg),
         paste("\nPercentage or enhancers mutated in at least 1 donor: ", perc)
       ),
       x = "", y = "number of donors")+
  theme_light()+
  theme(axis.text.x = element_blank())+
  scale_y_continuous(limits = c(0,30))
# print(b1)
#ggsave(paste(OUT_FOLDER_PLOTS, paste("Stat_enhancers", marker, "all_enhancers", mut, sep = "_"), ".png", sep = ""), plot = b1, device = "png", width = 10, height = 7)


# P-values
print("Wilcox - n_donors")
print(wilcox.test(stat_enh[,1], stat_enh[,2])[["p.value"]])
print(wilcox.test(stat_enh[,1], stat_enh[,3])[["p.value"]])

print("T-test - n_donors")
print(t.test(stat_enh[,1], stat_enh[,2])[["p.value"]])
print(t.test(stat_enh[,1], stat_enh[,3])[["p.value"]])

avg <- str_flatten(paste(names(apply(stat_don, 2, mean)), ":", round(apply(stat_don, 2, mean),2), sep = ""), 
                   collapse = "  ")
perc <- apply(stat_don, 2, function(x) round(sum(x != 0) / length(unique(SSMs$icgc_donor_id)) * 100))
perc <- str_flatten(paste(names(perc), ":", perc, " %", sep=""), collapse = "  ")
b2 <- stat_don %>%
  pivot_longer(everything(), names_to = "cluster", values_to = "value") %>%
  ggplot(., aes(x = cluster, y = log10(value), fill = cluster))+
  geom_boxplot()+
  labs(title = paste("Number of enahncers hit in each donor", sep=""), 
       subtitle = paste(
         paste("Average number of enhancers hit in a donor: ", avg), 
         paste("\nPercentage of donors with at least one enhancer mutated : ", perc)),
       x = "", y = "log10(number of enhancers hit)")+
  theme_light()+
  theme(axis.text.x = element_blank())
#print(b2)
#ggsave(paste(OUT_FOLDER_PLOTS, paste("Stat_donors", marker, "all_enhancers", mut, sep = "_"), ".png", sep = ""), plot = b1, device = "png", width = 10, height = 7)

# P-values
print(wilcox.test(stat_don[,1], stat_don[,2])[["p.value"]])
print(wilcox.test(stat_don[,1], stat_don[,3])[["p.value"]])
}


##

### Extracting SNVs coordinates for RegulomeDB ###

# RegulomeDB has 0-based coordinates
# ICGC has 1-based coordinates - https://docs.icgc.org/submission/guide/icgc-simple-somatic-mutation-format/

enh_SSMs_ctip <- read_tsv(paste0(OUT_FOLDER_DATA, "Table_enh_SSMs_CtIP.all_overlaps.", WIN, "kb_WIN.tsv"))
enh_SSMs_grhl <- read_tsv(paste0(OUT_FOLDER_DATA, "Table_enh_SSMs_GRHL.all_overlaps.", WIN, "kb_WIN.tsv"))


snvs_coords_ctip <- enh_SSMs_ctip %>% dplyr::select(seqnames_sbj, start_sbj, end_sbj) %>%
  filter(!duplicated(.)) %>%
  mutate(length = abs(end_sbj - start_sbj)) 
# Retain only SNVs (length = 0)
snvs_coords_ctip <- snvs_coords_ctip[snvs_coords_ctip$length == 0, ] 
# convert to 0-based system
snvs_coords_ctip$start_sbj <- snvs_coords_ctip$start_sbj - 1
snvs_coords_ctip$length <- NULL


snvs_coords_grhl <- enh_SSMs_grhl %>% dplyr::select(seqnames_sbj, start_sbj, end_sbj) %>%
  filter(!duplicated(.)) %>%
  mutate(length = abs(end_sbj - start_sbj)) 

# Retain only SNVs (length = 0)
snvs_coords_grhl <- snvs_coords_grhl[snvs_coords_grhl$length == 0, ] 

# convert to 0-based system
snvs_coords_grhl$start_sbj <- snvs_coords_grhl$start_sbj - 1
snvs_coords_grhl$length <- NULL


# Save coords in chr:start-end for RegulomeDB
str_c(str_c(snvs_coords_ctip$seqnames_sbj, snvs_coords_ctip$start_sbj, sep = ":"), snvs_coords_ctip$end_sbj, sep = "-") %>% 
  data.frame(.) #%>% write_tsv(., paste(OUT_FOLDER_DATA, "SNVs_coords_0_based.overlapping_with_enh_CtIP.for_regulomeDB.tsv", sep = ""), col_names = F)
str_c(str_c(snvs_coords_grhl$seqnames_sbj, snvs_coords_grhl$start_sbj, sep = ":"), snvs_coords_grhl$end_sbj, sep = "-") %>% 
  data.frame(.) #%>% write_tsv(., paste(OUT_FOLDER_DATA, "SNVs_coords_0_based.overlapping_with_enh_GRHL.for_regulomeDB.tsv", sep = ""), col_names = F)



##


## STSMs

# Position of Start/End bkpts within regions
df_var_stats_all <- list("CtIP" = NA, "GRHL" = NA)


for(marker in names(enh_all)){
  
  ### Pre/process input data
  marker_enh <- enh_all[[marker]] 
  
  ### Extend regions of 'WIN' from summit
  marker_enh$start <- marker_enh$summit - WIN
  marker_enh$end <- marker_enh$summit + WIN
  
  enh_gr <- makeGRangesFromDataFrame(marker_enh, keep.extra.columns = T)
  
  
  ### Find SNVs / Enhancers overlaps 
  hits_obj_enh <- findOverlaps(query=enh_gr, subject=STSMs_gr)
  STSMs_sbj <- data.frame(STSMs_gr[subjectHits(hits_obj_enh)])
  colnames(STSMs_sbj)[1:5] <- paste(colnames(STSMs_sbj)[1:5], "sbj", sep = "_")
  enh_SSMs <- cbind(data.frame(enh_gr[queryHits(hits_obj_enh)]), STSMs_sbj) 
  
  # By enhancer
  df_var_stats <- data.frame(name = unique(enh_SSMs$name), 
                             n_bkpts = NA, 
                             n_single_bkpts = NA,
                             n_within_stsms = NA)
  
  for(i in 1:length(unique(enh_SSMs$name))){
    # Variants overlapping each enhancer - ACROSS DONORS
    enh <- unique(enh_SSMs$name)[i]
    SSMs_sub <- enh_SSMs[enh_SSMs$name == enh, ]
    
    df_var_stats[df_var_stats$name == enh, "n_bkpts"] <- dim(SSMs_sub)[1]
    
    # For each enhancer, are there SVs with 2 bkpts within the enhancer?
    n_unique_vars <- table(SSMs_sub$sv_id)
    
    # For each var, only 1 instance. Only start or end bkpt falls within enhancer
    if(any(n_unique_vars > 1) == F){
      df_var_stats[df_var_stats$name == enh, "n_single_bkpts"] <- dim(SSMs_sub)[1]
      df_var_stats[df_var_stats$name == enh, "n_within_stsms"] <- 0
    } else if(any(n_unique_vars > 1) == T){
      df_var_stats[df_var_stats$name == enh, "n_single_bkpts"] <- sum(n_unique_vars == 1)
      
      n_single <- 0
      n_within <- 0
      
      # Which stsms have more than 1 bkpt
      SSMs_sub <- SSMs_sub %>% filter(sv_id %in% names(n_unique_vars[n_unique_vars > 1]))
      non_unique_vars <- unique(SSMs_sub$sv_id)
      for(i in 1:length(non_unique_vars)){
        one_var_df <- SSMs_sub[SSMs_sub$sv_id == non_unique_vars[i], ]
        # How many unique sv_id_x_class (from/to) x icgc_donor_id 
        df_stat <- one_var_df %>% group_by(icgc_donor_id) %>% summarise(n_sv_id_x_class = length(unique(sv_id_x_class)))
        n_within <- n_within + sum(df_stat$n_sv_id_x_class == 2)
        n_single <- n_single + sum(df_stat$n_sv_id_x_class == 1)
      }
      
      df_var_stats[df_var_stats$name == enh, "n_within_stsms"] <- n_within
      df_var_stats[df_var_stats$name == enh, "n_single_bkpts"] <- df_var_stats[df_var_stats$name == enh, "n_single_bkpts"]+n_single
    }
  }
  
  df_var_stats_all[[marker]] <- df_var_stats
}


for(marker in names(enh_all)){
  var_stats_marker <- df_var_stats_all[[marker]]
  print(apply(var_stats_marker[,-1], MARGIN = 2, FUN = sum))
  
  print(sprintf("Percentage of single bkpts: %.0f%%", 
                apply(var_stats_marker[,-1], MARGIN = 2, FUN = sum)[2] / apply(var_stats_marker[,-1], MARGIN = 2, FUN = sum)[1]*100))
  print(sprintf("Percentage of within bkpts: %.0f%%", 
                apply(var_stats_marker[,-1], MARGIN = 2, FUN = sum)[3] / apply(var_stats_marker[,-1], MARGIN = 2, FUN = sum)[1]*100))
}

aa <- left_join(var_stats_marker[var_stats_marker$n_within_stsms != 0, ], enh_SSMs[, c("name", "cluster")], by = "name")
aa <- aa[!duplicated(aa), ]
table(aa$cluster)



# Deletions spanning the full enhancer
deletions <- STSMs %>% filter(variant_type == "deletion")
length(unique(deletions$sv_id))
dim(deletions)[1]

deletions_from <- deletions %>% filter(class == "from") %>%
  select(-c(chrom_bkpt, end, range, class))
deletions_to <- deletions %>% filter(class == "to") %>%
  select(sv_id, start)
deletions_coords <- left_join(deletions_from, deletions_to, by = "sv_id")
colnames(deletions_coords)[12:13] <- c("start", "end")

# Length of deletions
data.frame(x = deletions_coords$end - deletions_coords$start) %>%
  ggplot(., aes(x = x))+
  geom_histogram()+
  scale_x_continuous(limits = c(0, 5.0e+07))
summary(deletions_coords$end - deletions_coords$start)

deletions_gr <- makeGRangesFromDataFrame(deletions_coords, keep.extra.columns = T)


for(marker in names(enh_all)){
  marker_enh <- enh_all[[marker]] 
  marker_enh$start <- marker_enh$summit - WIN
  marker_enh$end <- marker_enh$summit + WIN
  
  enh_gr <- makeGRangesFromDataFrame(marker_enh, keep.extra.columns = T)
  
  # if type is within, the query interval must be wholly contained within the subject interval. 
  hits_obj <- findOverlaps(query = enh_gr, subject = deletions_gr, 
                           type = "within", minoverlap = 1)
  del_sbj <- data.frame(deletions_gr[subjectHits(hits_obj)])
  colnames(del_sbj)[1:5] <- paste(colnames(del_sbj)[1:5], "sbj", sep = "_")
  enh_dels <- cbind(data.frame(enh_gr[queryHits(hits_obj)]), del_sbj) 
  
  sum(enh_dels$end_sbj - enh_dels$end < 0) # should be 0 - STSM end always after enhancer end
  sum(enh_dels$start_sbj - enh_dels$start > 0) # should be 0 - STSM start always before enhancer start
  
  cat("\n")
  print(marker)
  print(sprintf("Number of deletions spanning a full enhancer: %.0f - %.1f%%", length(unique(enh_dels$sv_id)), 
                length(unique(enh_dels$sv_id))/dim(deletions_coords)[1]*100))  
  print(sprintf("Number of enhancers fully deleted: %.0f - %.1f%%", length(unique(enh_dels$name)), 
                length(unique(enh_dels$name))/dim(marker_enh)[1]*100))
  print(sprintf("Number of donors hit: %.0f",length(unique(enh_dels$icgc_donor_id))))
  
  # Number of enhancers removed by each deletion
  n_del_enh <- enh_dels %>% group_by(sv_id) %>% 
    dplyr::summarise(n_enh_deleted = length(unique(name))) %>%
    arrange(-n_enh_deleted)
  
  print("Number of enhancers removed by each deletion:")
  print(summary(n_del_enh$n_enh_deleted))
  print("Number of deletions removing only one enhancer:")
  print(dim(n_del_enh[n_del_enh$n_enh_deleted == 1, ])[1])
  
  
  # Number of donors harboring an 'enhancer-specific' deletion
  enh_spec_del <- n_del_enh[n_del_enh$n_enh_deleted == 1, ]$sv_id
  n_enh_spec_dels <- enh_dels %>% filter(sv_id %in% enh_spec_del) %>% 
    group_by(icgc_donor_id) %>% 
    dplyr::summarise(n_enh_spec_dels = length(unique(sv_id)))
  n_enh_spec_dels <- rbind(n_enh_spec_dels,
                           data.frame("icgc_donor_id" = unique(STSMs$icgc_donor_id[!STSMs$icgc_donor_id %in% n_enh_spec_dels$icgc_donor_id]),
                                      "n_enh_spec_dels" = 0))
  
  print("Number of 'enhancer-specific' dels x donor:")
  print(summary(n_enh_spec_dels$n_enh_spec_dels))
  
}



