
library(tidyverse)

STSMs_ctip <- read_tsv("/Users/ieo6983/Desktop/enhancers_project_old/per_ivan/STSMs_overlapping_CTIP_enhancers.tsv")
STSMs_grhl <- read_tsv("/Users/ieo6983/Desktop/enhancers_project_old/per_ivan/STSMs_overlapping_GRHL_enhancers.tsv")

STSMs_all <- read_tsv("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/genomics/pre_processed_ICGC/structural_somatic_mutation.preprocessed.tsv")
STSMs_raw <- read_tsv("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/genomics/raw_ICGC/structural_somatic_mutation.tsv")


##


STSMs_all_to_merge <- STSMs_all %>% dplyr::select(., chrom, chrom_bkpt, icgc_donor_id, variant_type, sv_id, class)
colnames(STSMs_all_to_merge)[1:2] <- c("seqnames_sbj2", "start_sbj2")
STSMs_all_to_merge <- STSMs_all_to_merge %>% mutate(., end_sbj2 = .$start_sbj2, .after = start_sbj2)
STSMs_all_to_merge$seqnames_sbj2 <- paste0("chr", STSMs_all_to_merge$seqnames_sbj2)

# Are all sv_id present?
sum(!STSMs_ctip$sv_id %in% STSMs_all_to_merge$sv_id)
sum(!STSMs_grhl$sv_id %in% STSMs_all_to_merge$sv_id)


##


## Add second bkpt - in 2 steps

# get unique bkpts from overlaps df
STSMs_ctip1 <- STSMs_ctip[!duplicated(STSMs_ctip$sv_id), ]
STSMs_ctip2 <- STSMs_ctip[duplicated(STSMs_ctip$sv_id), ]


ctip_to_merge <- STSMs_all_to_merge %>% 
  semi_join(., STSMs_ctip, by = "sv_id") %>%
  anti_join(., STSMs_ctip, by = c("sv_id", "class"))

left_join()


##

#CtIP
STSMs_ctip_new <- data.frame()

for(i in 1:dim(STSMs_ctip)[1]){
  row_1 <- STSMs_ctip[i, ]
  row_2 <- STSMs_all_to_merge[ (STSMs_all_to_merge$sv_id == row_1$sv_id & !STSMs_all_to_merge$class == row_1$class) , ]
  connected_rows <- cbind(row_1, row_2)
  
  STSMs_ctip_new <- rbind(STSMs_ctip_new, connected_rows)
}
colnames(STSMs_ctip_new) <- c(
  colnames(STSMs_ctip_new)[1:16], 
  c("icgc_donor_id2", "variant_type2", "sv_id2", "class2")
)

sum(!STSMs_ctip_new$sv_id == STSMs_ctip_new$sv_id2)
sum(STSMs_ctip_new$class == STSMs_ctip_new$class2)


#GRHL2
STSMs_grhl_new <- data.frame()

for(i in 1:dim(STSMs_grhl)[1]){
  row_1 <- STSMs_grhl[i, ]
  row_2 <- STSMs_all_to_merge[ (STSMs_all_to_merge$sv_id == row_1$sv_id & !STSMs_all_to_merge$class == row_1$class) , ]
  connected_rows <- cbind(row_1, row_2)
  
  STSMs_grhl_new <- rbind(STSMs_grhl_new, connected_rows)
}
colnames(STSMs_grhl_new) <- c(
  colnames(STSMs_grhl_new)[1:16], 
  c("icgc_donor_id2", "variant_type2", "sv_id2", "class2")
)

sum(!STSMs_grhl_new$sv_id == STSMs_grhl_new$sv_id2)
sum(STSMs_grhl_new$class == STSMs_grhl_new$class2)


##


#write_tsv(STSMs_ctip_new, "/Users/ieo6983/Desktop/enhancers_project_old/per_ivan/STSMs_overlapping_CTIP_enhancers.both_bkpts.tsv")
#write_tsv(STSMs_grhl_new, "/Users/ieo6983/Desktop/enhancers_project_old/per_ivan/STSMs_overlapping_GRHL_enhancers.both_bkpts.tsv")
