
#'Annotate SVs as damaging or non-damaging
#'
#'@Deletion is damaging if 1 bkpt falls within E - P region and the other one falls outside. 
#'Deletions in which both bpts fall within E - P region bring enhancer and promoter closer.
#'@Interchromosomal-translocation is damaging if 1 bkpt falls within E - P region and the other one falls outside.
#'two bkpts of the same SV cannot fall within the same loops, because they fall in different chromosomes.
#'@Inversion is damaging is damaging if 1 bkpt falls within E - P region and the other one falls outside
#'and the distance 'outer' bkpt - bin is (1) greater than the distance bin - 'inner' bkpt, (2) greater than distance (E-P) / 2 

library(tidyverse)
library(GenomicRanges)
library(assertthat)

SEED <- 4321
set.seed(SEED)

#MARKERS <- c("CtIP", "GRHL")
marker <- "GRHL"
WIN <- 3000 # enhancers window that was used to compute overlap with SVs
loops_kb <- 2
naive <- F
save_SVs_anno <- F

path_SVs <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/genomics/raw_ICGC/structural_somatic_mutation.tsv")
path_loops <- fs::path(paste0("/Users/ieo6983/Desktop/fragile_enhancer_clinical/results/integrated/OLD/", loops_kb, "kb/data/2kb_Unified_table.SCR_plus_KD_counts.all_anno_loops.ENH_DEGs_any.tsv"))
path_enhancers <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/others/Cluster_GRHL_Enh_All.txt")
path_results_data <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/results/ICGC/OLD/damaging_variants_annotation/SVs/")  

                              
##

# TODO @OPTIMIZE: naive function. Could make it handle more cases. 
anno_distinct_loops <- function(df, col1, col2){
  if(is.string(col1) != T | is.string(col2) != T){
    print("Params must be strings")
  }

  col1_str <- df[, col1]  
  col2_str <- df[, col2]
  
  # Given that SVs are right-oriented, loop 2 comes AFTER loop 1 (at least for 1 bin)
  col1_split <- str_split(col1_str, pattern = "_", simplify = T)
  col2_split <- str_split(col2_str, pattern = "_", simplify = T)
  
  # loop2-start > loop1-end
  cond1 <- (as.numeric(col2_split[, 2])) > (as.numeric(col1_split[, 3]))
  
  # If NA | TRUE = damaging. If FALSE = false_damaging
  df$distinct_loops <- F
  df$distinct_loops[ (is.na(cond1) | cond1 == TRUE) ] <- T
  
  return(df)
}


##


# Read SVs, loops and enhancers
SVs_full <- read_tsv(path_SVs)
loops <- read_tsv(path_loops)
enh <- read_tsv(path_enhancers, comment = "#", col_names = c("chrom", "start", "end", "cluster"))
enh$summit <- enh$end
enh$name <- str_c(enh$chrom, enh$summit, sep=":")

print("Types of structural variants:")
print(table(SVs_full$variant_type))
# Types are: 
# "deletion", "inversion", "interchromosomal rearrangement - unknown type"
# "intrachromosomal rearrangement with inverted orientation", "intrachromosomal rearrangement with non-inverted orientation"
# "tandem duplication"

## Simplify variant types (leaving out variant_type = "tandem duplication")
variant_types_complex <- c("deletion", "inversion", "interchromosomal rearrangement - unknown type")
variant_types_simple <- c("deletion", "inversion", "interchromosomal-translocation")
types_conv <- data.frame("complex" = variant_types_complex, "simple" = variant_types_simple)

SVs_full <- SVs_full %>% left_join(., types_conv, by = c("variant_type" = "complex")) %>% 
  rename(., "simple" = "variant_type_simple") %>% relocate(., variant_type_simple, .after = variant_type)
SVs_full$chr_from <- paste0("chr", SVs_full$chr_from)
SVs_full$chr_to <- paste0("chr", SVs_full$chr_to)
# Reduce to useful columns
SVs <- SVs_full[, 1:22] %>% dplyr::select(., -c(placement, interpreted_annotation, chr_from_flanking_seq, chr_to_flanking_seq))

## Define E - P regions to be intersected with SVs
loops_enh_deg <- loops %>% dplyr::filter(., ( !is.na(name1) & !is.na(gene_name2) ) | 
                           ( !is.na(gene_name1) & !is.na(name2) ) )
# end1 is always < end2 & end1 is always < start2
loops_enh_deg$EP_chrom <- loops_enh_deg$seqnames1
loops_enh_deg$EP_start <- loops_enh_deg$start1
loops_enh_deg$EP_end <- loops_enh_deg$end2
loops_enh_deg$unique_loop_id <- paste(loops_enh_deg$EP_chrom, loops_enh_deg$EP_start, loops_enh_deg$EP_end, sep = "_")
loops_enh_deg_gr <- loops_enh_deg %>% makeGRangesFromDataFrame(., keep.extra.columns = T, 
                                           seqnames.field = "EP_chrom", start.field = "EP_start", end.field = "EP_end")

## Annotate SVs
print("All types: ") 
print(unique(SVs$variant_type_simple))

SVs_annotated <- data.frame()
for(type in unique(SVs$variant_type_simple)){
  
  # tandem duplication - not considered
  if(is.na(type)){
    next
  }
  
  if(type == "deletion"){
    dels <- SVs %>% dplyr::filter(., variant_type_simple == type)
    print(paste0("--- Annotating variants of type: ", type, " ---"))
    print(paste0("Total number of ", type, " : ", dim(dels)[1]))

    if(naive == T){
      # Check whether bkpt 1 and bkpt 2 of overlap with one E - P region
      dels_from <- makeGRangesFromDataFrame(dels[, c("chr_from", "chr_from_bkpt")], seqnames.field = "chr_from", start.field = "chr_from_bkpt", end.field = "chr_from_bkpt")
      dels_to <- makeGRangesFromDataFrame(dels[, c("chr_to", "chr_to_bkpt")], seqnames.field = "chr_to", start.field = "chr_to_bkpt", end.field = "chr_to_bkpt")
      
      filt_from <- !is.na(findOverlaps(dels_from, loops_enh_deg_gr, select = "arbitrary"))
      filt_to <- !is.na(findOverlaps(dels_to, loops_enh_deg_gr, select = "arbitrary"))
      
      dels$ovrlp_from <- F
      dels$ovrlp_to <- F
      dels[filt_from, ]$ovrlp_from <- T
      dels[filt_to, ]$ovrlp_to <- T
      
      # Is deletion damaging? How many damaging deletions over total?
      dam_vec <- dels$ovrlp_from + dels$ovrlp_to
      dels$damaging <- F
      dels[dam_vec == 1, ]$damaging <- T
      
      print(paste0("Damaging deletions: ", sum(dels$damaging), " over ", length(dels$damaging), 
                   " (", round(sum(dels$damaging) / length(dels$damaging),2)*100, "%)")) 
      
      message(paste0("Potentially, ", sum(dam_vec == 2), " false negatives."))
      cat("\n") 
    }
      
    ##
    
    # Annotate SVs to loops and damaging effect 
    dels <- SVs %>% dplyr::filter(., variant_type_simple == type)
    dels_from <- makeGRangesFromDataFrame(dels, seqnames.field = "chr_from", start.field = "chr_from_bkpt", end.field = "chr_from_bkpt", keep.extra.columns = T)
    dels_to <- makeGRangesFromDataFrame(dels, seqnames.field = "chr_to", start.field = "chr_to_bkpt", end.field = "chr_to_bkpt", keep.extra.columns = T)
    df_from <- dels_from; df_to <- dels_to;
    
    hits_from <- findOverlaps(query=df_from, subject=loops_enh_deg_gr)
    overlaps_from <- cbind(data.frame(df_from[queryHits(hits_from)]),
                           data.frame("unique_loop_id_from" = mcols(loops_enh_deg_gr[subjectHits(hits_from)])[, "unique_loop_id"]))
    overlaps_from <- overlaps_from %>% dplyr::rename(.,  "chr_from" = "seqnames", "chr_from_bkpt" = "start") %>% 
      dplyr::select(., -c(end, width, strand)) %>% dplyr::relocate(., c(chr_from, chr_from_bkpt), .before = chr_from_strand)    
    
    hits_to <- findOverlaps(query=dels_to, subject=loops_enh_deg_gr)
    overlaps_to <- cbind(data.frame(dels_to[queryHits(hits_to)]),
                         data.frame("unique_loop_id_to" = mcols(loops_enh_deg_gr[subjectHits(hits_to)])[, "unique_loop_id"]))
    overlaps_to <- overlaps_to %>% dplyr::rename(.,  "chr_to" = "seqnames", "chr_to_bkpt" = "start") %>% 
      dplyr::select(., -c(end, width, strand)) %>% dplyr::relocate(., c(chr_to, chr_to_bkpt), .before = chr_to_strand)    
    
    # full_join joins by all columns in common (only unique_loop_id_from/to is not in common)
    overlaps_join <- full_join(overlaps_from, overlaps_to, relationship = "many-to-many") %>% suppressMessages()
    
    # Only SVs overlapping: (1) one loop - NA, or (2) two different loops are damaging   
    overlaps_join <- anno_distinct_loops(df=overlaps_join, col1="unique_loop_id_from", col2="unique_loop_id_to") 
    
    # For each deletion, compute in how many cases it was damaging 
    dam_vec <- overlaps_join %>% group_by(., sv_id) %>% dplyr::summarise(., "how_many_dam" = sum(distinct_loops)) %>% .[,"how_many_dam"]
    print(paste0("Damaging deletions: ", sum(dam_vec >= 1), " over ", length(unique(dels$sv_id)), 
                 " (", round(sum(dam_vec >= 1) / length(unique(dels$sv_id)),2)*100, "%)")) 
    cat("\n") 
    
    # For which loop is SV damaging? (from or to)
    overlaps_join$dam_from <- ( !is.na(overlaps_join$unique_loop_id_from) ) & ( overlaps_join$distinct_loops == T)
    overlaps_join$dam_to <- ( !is.na(overlaps_join$unique_loop_id_to) ) & ( overlaps_join$distinct_loops == T)
    overlaps_join <- overlaps_join %>% rename(., "distinct_loops" = "damaging")

    # Save annotated SVs
    SVs_annotated <- rbind(SVs_annotated, overlaps_join)
  } 
  
  else if(type == "interchromosomal-translocation"){
    trans <- SVs %>% dplyr::filter(., variant_type_simple == type)
    print(paste0("--- Annotating variants of type: ", type, " ---"))
    print(paste0("Total number of ", type, ": ", dim(trans)[1])) 
   
    if(naive == T){
      # Check whether bkpt 1 and bkpt 2 of overlap with one E - P region
      trans_from <- makeGRangesFromDataFrame(trans[, c("chr_from", "chr_from_bkpt")], seqnames.field = "chr_from", start.field = "chr_from_bkpt", end.field = "chr_from_bkpt")
      trans_to <- makeGRangesFromDataFrame(trans[, c("chr_to", "chr_to_bkpt")], seqnames.field = "chr_to", start.field = "chr_to_bkpt", end.field = "chr_to_bkpt")
      
      filt_from <- !is.na(findOverlaps(trans_from, loops_enh_deg_gr, select = "arbitrary")) %>% suppressWarnings()
      filt_to <- !is.na(findOverlaps(trans_to, loops_enh_deg_gr, select = "arbitrary")) %>% suppressWarnings()
      
      trans$ovrlp_from <- F
      trans$ovrlp_to <- F
      trans[filt_from, ]$ovrlp_from <- T
      trans[filt_to, ]$ovrlp_to <- T
      
      # Is translocation damaging? How many damaging deletions over total?
      dam_vec <- trans$ovrlp_from + trans$ovrlp_to
      trans$damaging <- F
      # Both 1 and 2 denote a damaging translocation. This because the two bkpts cannot fall within the same loop (different chroms)
      # Hence if dam_vec == 2, the 2 bkpts fall within 2 different loops
      trans[!dam_vec == 0 , ]$damaging <- T
      
      print(paste0("Damaging translocations: ", sum(trans$damaging), " over ", length(trans$damaging), 
                   " (", round(sum(trans$damaging) / length(trans$damaging),2)*100, "%)")) 
      cat("\n") 
    }
    
    ##
    
    # Annotate SVs to loops and damaging effect 
    trans <- SVs %>% dplyr::filter(., variant_type_simple == type)
    trans_from <- makeGRangesFromDataFrame(trans, seqnames.field = "chr_from", start.field = "chr_from_bkpt", end.field = "chr_from_bkpt", keep.extra.columns = T)
    trans_to <- makeGRangesFromDataFrame(trans, seqnames.field = "chr_to", start.field = "chr_to_bkpt", end.field = "chr_to_bkpt", keep.extra.columns = T)
    df_from <- trans_from; df_to <- trans_to;
    
    hits_from <- findOverlaps(query=df_from, subject=loops_enh_deg_gr)
    overlaps_from <- cbind(data.frame(df_from[queryHits(hits_from)]),
                           data.frame("unique_loop_id_from" = mcols(loops_enh_deg_gr[subjectHits(hits_from)])[, "unique_loop_id"]))
    overlaps_from <- overlaps_from %>% dplyr::rename(.,  "chr_from" = "seqnames", "chr_from_bkpt" = "start") %>% 
      dplyr::select(., -c(end, width, strand)) %>% dplyr::relocate(., c(chr_from, chr_from_bkpt), .before = chr_from_strand)    
    
    hits_to <- findOverlaps(query=df_to, subject=loops_enh_deg_gr) %>% suppressWarnings() 
    overlaps_to <- cbind(data.frame(df_to[queryHits(hits_to)]),
                         data.frame("unique_loop_id_to" = mcols(loops_enh_deg_gr[subjectHits(hits_to)])[, "unique_loop_id"]))
    overlaps_to <- overlaps_to %>% dplyr::rename(.,  "chr_to" = "seqnames", "chr_to_bkpt" = "start") %>% 
      dplyr::select(., -c(end, width, strand)) %>% dplyr::relocate(., c(chr_to, chr_to_bkpt), .before = chr_to_strand) 
    
    # full_join joins by all columns in common (only unique_loop_id_from/to is not in common)
    overlaps_join <- full_join(overlaps_from, overlaps_to, relationship = "many-to-many") %>% suppressMessages()
    
    # All translocations within overlaps_join should be damaging
    dam_vec <- is.na(overlaps_join$unique_loop_id_from) & is.na(overlaps_join$unique_loop_id_to)
    if(sum(dam_vec) == 0){
      print(paste0("Damaging interchromosomal-translocations: ", length(unique(overlaps_join$sv_id)), " over ", length(unique(trans$sv_id)), 
                   " (", round( length(unique(overlaps_join$sv_id)) / length(unique(trans$sv_id)), 2)*100, "%)"))
      cat("\n")
      
      overlaps_join$damaging <- TRUE
      # For which loop is the SV damaging? (from or to)
      overlaps_join$dam_from <- !is.na(overlaps_join$unique_loop_id_from)
      overlaps_join$dam_to <- !is.na(overlaps_join$unique_loop_id_to)
  
      # Save annotated SVs 
      SVs_annotated <- rbind(SVs_annotated, overlaps_join)
    } else { message(paste0("Check ", type)) }    
  }
  
  else if(type == "inversion"){
    invs <- SVs %>% dplyr::filter(., variant_type_simple == type)
    print(paste0("--- Annotating variants of type: ", type, " ---"))
    print(paste0("Total number of ", type, ": ", dim(invs)[1])) 
    
    # Check whether bkpt 1 and bkpt 2 of overlap with one E - P region
    invs_from <- makeGRangesFromDataFrame(invs, seqnames.field = "chr_from", start.field = "chr_from_bkpt", end.field = "chr_from_bkpt", keep.extra.columns = T)
    invs_to <- makeGRangesFromDataFrame(invs, seqnames.field = "chr_to", start.field = "chr_to_bkpt", end.field = "chr_to_bkpt", keep.extra.columns = T)
    df_from <- invs_from; df_to <- invs_to;
    
    # Overlap separately from and to, keep which is the loop_id, then merge the 2 datasets keeping this info for both SVs coords 
    hits_from <- findOverlaps(query=df_from, subject=loops_enh_deg_gr)
  
    overlaps_from <- cbind(data.frame(df_from[queryHits(hits_from)]),
                           data.frame("unique_loop_id_from" = mcols(loops_enh_deg_gr[subjectHits(hits_from)])[, "unique_loop_id"]))
    overlaps_from <- overlaps_from %>% dplyr::rename(.,  "chr_from" = "seqnames", "chr_from_bkpt" = "start") %>% 
      dplyr::select(., -c(end, width, strand)) %>% dplyr::relocate(., c(chr_from, chr_from_bkpt), .before = chr_from_strand)    
    
    hits_to <- findOverlaps(query=df_to, subject=loops_enh_deg_gr) %>% suppressWarnings() 
    overlaps_to <- cbind(data.frame(df_to[queryHits(hits_to)]),
                         data.frame("unique_loop_id_to" = mcols(loops_enh_deg_gr[subjectHits(hits_to)])[, "unique_loop_id"]))
    overlaps_to <- overlaps_to %>% dplyr::rename(.,  "chr_to" = "seqnames", "chr_to_bkpt" = "start") %>% 
      dplyr::select(., -c(end, width, strand)) %>% dplyr::relocate(., c(chr_to, chr_to_bkpt), .before = chr_to_strand) 
    
    # full_join joins by all columns in common (only unique_loop_id_from/to is not in common)
    overlaps_join <- full_join(overlaps_from, overlaps_to, relationship = "many-to-many") %>% suppressMessages()
    
    ##
    
    # Handle options 
    
    #' @Option 1: none of from/to SV coords overlap with a loop
    #' These SVs are not present in overlaps_join, which contains only SVs-loops matches
    
    #' @Option 2: one of the two SV coords (either from or to) overlaps with a loop
    sum(!invs$chr_from_bkpt <= invs$chr_to_bkpt) # always T - 2nd bkpt always comes after 1st bkpt
    opt2_vec<- (!is.na(overlaps_join$unique_loop_id_from)) + (!is.na(overlaps_join$unique_loop_id_to))
    opt_2 <- overlaps_join %>% dplyr::filter(., 
                                             (!is.na(unique_loop_id_from) & is.na(unique_loop_id_to)) | (is.na(unique_loop_id_from) & !is.na(unique_loop_id_to)) )
    
    # 1st bkpt overlaps with loop
    # Compute distances: inner bkpt_from - loop_bin2 & loop_bin2 - outer bkpt_to
    opt2_from <- opt_2[!is.na(opt_2$unique_loop_id_from), ] %>% 
      left_join(., loops_enh_deg, by = c("unique_loop_id_from" = "unique_loop_id")) %>% 
      mutate(., "inner_dist" = ( (start2+(loops_kb*1000/2)) - chr_from_bkpt) ) %>%
      mutate(., "outer_dist" = ( chr_to_bkpt -  (start2+(loops_kb*1000/2))) ) %>%
      mutate(., "outer_minus_inner" = outer_dist - inner_dist) %>%
      mutate(., "damaging" = ifelse(outer_minus_inner > ((EP_end - EP_start - 2000) / 2), TRUE, FALSE)) %>% 
      dplyr::select(., -c(seqnames1:outer_minus_inner))
    
    # Compute distances: inner bkpt_to - loop_bin1 & loop_bin1 - outer bkpt_from
    opt2_to <- opt_2[!is.na(opt_2$unique_loop_id_to), ] %>% 
        left_join(., loops_enh_deg, by = c("unique_loop_id_to" = "unique_loop_id")) %>%  
        mutate(., "inner_dist" = ( (chr_to_bkpt - (start1+(loops_kb*1000/2))) ) ) %>% 
        mutate(., "outer_dist" = ( (start1+(loops_kb*1000/2)) - chr_from_bkpt) ) %>%
        mutate(., "outer_minus_inner" = outer_dist - inner_dist) %>%
        mutate(., "damaging" = ifelse(outer_minus_inner > ((EP_end - EP_start - 2000) / 2), TRUE, FALSE))%>% 
      dplyr::select(., -c(seqnames1:outer_minus_inner))
    
    opt_2_new <- rbind(opt2_from, opt2_to)
    opt_2_new$dam_from <- ( !is.na(opt_2_new$unique_loop_id_from) & (opt_2_new$damaging == T) )
    opt_2_new$dam_to <- ( !is.na(opt_2_new$unique_loop_id_to) & (opt_2_new$damaging == T) )
    
    #' @Option 3: both SV coords (from and to) overlap with a loop
    opt_3 <- overlaps_join %>% dplyr::filter(., (!is.na(unique_loop_id_from) & !is.na(unique_loop_id_to) ) ) 
    opt_3 <- anno_distinct_loops(df=opt_3, col1 = "unique_loop_id_from", "unique_loop_id_to")
    opt_3$dam_from <- F
    opt_3$dam_to <- F
    
    # Check whether inversion is damaging for one of the 2 loops
    opt_3_not <- opt_3 %>% dplyr::filter(., distinct_loops == F) 
    opt_3_dam <- opt_3 %>% dplyr::filter(., distinct_loops == T)
    
    # Compute distances: inner bkpt_from - loop_bin2 & loop_bin2 - outer bkpt_to
    opt_3_dam <- opt_3_dam %>% 
      left_join(., loops_enh_deg, by = c("unique_loop_id_from" = "unique_loop_id")) %>% 
      mutate(., "inner_dist" = ( (start2+(loops_kb*1000/2)) - chr_from_bkpt) ) %>%
      mutate(., "outer_dist" = ( chr_to_bkpt -  (start2+(loops_kb*1000/2))) ) %>%
      mutate(., "outer_minus_inner" = outer_dist - inner_dist) %>%
      mutate(., "dam_from" = ifelse(outer_minus_inner > ((EP_end - EP_start - 2000) / 2), TRUE, FALSE)) %>%
      dplyr::select(., -(seqnames1:outer_minus_inner))
    # Compute distances: inner bkpt_to - loop_bin1 & loop_bin1 - outer bkpt_from
    opt_3_dam <- opt_3_dam %>% 
      left_join(., loops_enh_deg, by = c("unique_loop_id_to" = "unique_loop_id")) %>%  
      mutate(., "inner_dist" = ( (chr_to_bkpt - (start1+(loops_kb*1000/2))) ) ) %>% 
      mutate(., "outer_dist" = ( (start1+(loops_kb*1000/2)) - chr_from_bkpt) ) %>%
      mutate(., "outer_minus_inner" = outer_dist - inner_dist) %>%
      mutate(., "dam_to" = ifelse(outer_minus_inner > ((EP_end - EP_start - 2000) / 2), TRUE, FALSE)) %>%
      dplyr::select(., -(seqnames1:outer_minus_inner))
    
    opt_3_full <- rbind(opt_3_dam, opt_3_not)
    
    opt_3_new <- opt_3_full
    opt_3_new$damaging <- F
    opt_3_new[opt_3_new$dam_from + opt_3_new$dam_to != 0, ]$damaging <- T
    # Keep info: for which loop is the SV damaging? (from or to)
    opt_3_new <- opt_3_new %>% dplyr::select(, -c(distinct_loops)) 
      
    ##
    
    # Ri-merge everything 
    
    overlaps_join <- rbind(opt_2_new, opt_3_new) %>% 
      relocate(., damaging, .before = dam_from)
    print(paste0("Damaging inversions ", length(unique(overlaps_join[overlaps_join$damaging == T,]$sv_id)), " over ", length(unique(invs$sv_id)), 
                 " (", round( length(unique(overlaps_join[overlaps_join$damaging == T,]$sv_id)) / length(unique(invs$sv_id)), 2)*100, "%)"))
    cat("\n")
    
    # Save annotated SVs
    SVs_annotated <- rbind(SVs_annotated, overlaps_join)
  }
    
}


##


## Save output 

# SVs_annotated contains all SVs that overlapped with an E - P region, classified either as 'damaging' or 'non-damaging' 
# If damaging = TRUE, SV is damaging for both loops indicated (unique_loop_id_from and unique_loop_id_to)

if(save_SVs_anno == T){
  write_tsv(SVs_annotated, file = fs::path(path_results_data, paste0("Table_SVs_loops_GRHL.all_overlaps.", loops_kb, "kb_res.with_damaging_anno.tsv")))
}


##


#' @NOTES: 
#' 
#' @E-Pregions overlapping with SVs: how do we define them? 
#' For now, E-P region goes from the beginning of bin1 (start1) and the end of bin2 (end2).
#' However, enhancers and promoters were annotated to loops by considering +/- 1 bin. 
#' So maybe we should either 
#' - (1) extend E-P region to (start1 - 2kb) - (end2 + 2kb) or
#' - (2) define E-P regions by considering where the enhancer/promoter falls within the bin  
#' 
#' @Neglected-SVs: some SVs types were overlooked, namely:
#' "intrachromosomal rearrangement with inverted orientation", "intrachromosomal rearrangement with non-inverted orientation", "tandem duplication"
#' Should understand what the first 2 are and eventually include them.
#' 
#' @Inversions are annotated based on 'outer' and 'inner' distances between SVs bkpts and the _middle_ point of bin1 or bin2  
#' Ideally, we should compute the distance among bkpts and enhancers-summit and DEGs-TSSs
#' But this is complicated, because EACH loop can overlap BOTH 1 enhancer and 1 promoter 
#' Should use unbudled loops to handle this. In any case, the number of damagin mutations is pretty low.
#' 
#' @Outside, definition: 
#' For both deletions and inversions, we consider as damaging those SVs with a 1st bkpts 'inside' a loop and a 2nd 'outside'.
#' 'outside' can be a non-annotated genomic region, the same loop as the 1st bkpt, or _another loop_.
#' In the last scenario, the SV could still be damaging, depending on cases. 
#' For now, the SV is still damagin if 2nd bkpt falls in a _distinct_ loop (with bin1 and bin2 outside of the first E-P region)
#' There are other complex situations, such as when one bin falls within the first E-P region, and the other does not. 
#' These are not handled yet.  
#' 
#' @Annotation:
#' Now we are considering loops connecting enhancers to DEGs promoters. Expand to enhancer - any gene?




