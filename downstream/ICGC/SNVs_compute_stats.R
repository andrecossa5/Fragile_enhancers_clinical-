
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





