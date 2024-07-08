
#'Annotate SVs as damaging or non-damaging
#'
#'
#'@Deletion is damaging if 1 bkpt falls within E - P region and the other one falls outside. 
#'Deletions in which both bpts fall within E - P region bring enhancer and promoter closer.
#'@Interchromosomal-translocation is damaging if 1 bkpt falls within E - P region and the other one falls outside.
#'two bkpts of the same SV cannot fall within the same loops, because they fall in different chromosomes.

library(tidyverse)
library(GenomicRanges)
#library(assertthat)

SEED <- 4321
set.seed(SEED)

#MARKERS <- c("CtIP", "GRHL")
marker <- "GRHL"
WIN <- 3000 # enhancers window that was used to compute overlap with SVs
loops_kb <- 2
naive <- F

path_SVs <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/genomics/raw_ICGC/structural_somatic_mutation.tsv")
path_loops <- fs::path(paste0("/Users/ieo6983/Desktop/fragile_enhancer_clinical/results/integrated/", loops_kb, "kb/data/2kb_Unified_table.SCR_plus_KD_counts.all_anno_loops.ENH_DEGs_any.tsv"))
path_enhancers <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/others/Cluster_GRHL_Enh_All.txt")
  

##


# TODO: naive function. Could make it handle more cases. 
anno_true_damaging <- function(df, col1, col2){
  col1_str <- df[, col1]  
  col2_str <- df[, col2]
  
  # Given that SVs are right-oriented, loop 2 comes AFTER loop 1 (at least for 1 bin)
  col1_split <- str_split(col1_str, pattern = "_", simplify = T)
  col2_split <- str_split(col2_str, pattern = "_", simplify = T)
  
  # loop2-start > loop1-end
  cond1 <- (as.numeric(col2_split[, 2])) > (as.numeric(col1_split[, 3]))
  
  # If NA | TRUE = damaging. If FALSE = false_damaging
  df$true_dam <- F
  df$true_dam[ (is.na(cond1) | cond1 == TRUE) ] <- T
  
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
print(table(SVs$variant_type))
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
    overlaps_join <- anno_true_damaging(df=overlaps_join, col1="unique_loop_id_from", col2="unique_loop_id_to") 
    
    dam_vec <- overlaps_join %>% group_by(., sv_id) %>% dplyr::summarise(., "how_many_dam" = sum(true_dam)) %>% .[,"how_many_dam"]
    print(paste0("Damaging deletions: ", sum(dam_vec >= 1), " over ", length(unique(dels$sv_id)), 
                 " (", round(sum(dam_vec >= 1) / length(unique(dels$sv_id)),2)*100, "%)")) 
    cat("\n") 
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
      # Both 1 and 2 denote a damaging translocation. This becasue the two bkpts cannot fall within the same loop (different chroms)
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
    } else { message(paste0("Check ", type)) }    
  }
  
  else if(type == "inversion"){
    invs <- SVs %>% dplyr::filter(., variant_type_simple == type)
    print(paste0("--- Annotating variants of type: ", type, " ---"))
    print(paste0("Total number of ", type, ": ", dim(invs)[1])) 
    
    # Check whether bkpt 1 and bkpt 2 of overlap with one E - P region
    invs_from <- makeGRangesFromDataFrame(invs[, c("chr_from", "chr_from_bkpt")], seqnames.field = "chr_from", start.field = "chr_from_bkpt", end.field = "chr_from_bkpt")
    invs_to <- makeGRangesFromDataFrame(invs[, c("chr_to", "chr_to_bkpt")], seqnames.field = "chr_to", start.field = "chr_to_bkpt", end.field = "chr_to_bkpt")
    
    filt_from <- !is.na(findOverlaps(invs_from, loops_enh_deg_gr, select = "arbitrary")) %>% suppressWarnings()
    filt_to <- !is.na(findOverlaps(invs_to, loops_enh_deg_gr, select = "arbitrary")) %>% suppressWarnings()
    
    invs$ovrlp_from <- F
    invs$ovrlp_to <- F
    invs[filt_from, ]$ovrlp_from <- T
    invs[filt_to, ]$ovrlp_to <- T
    
    # Keeping overlapping loop information 
    # If I keep the info of which loop the SV is overlapping (unique_loop_id), I can re-connect the SVs with the loops table
    # From there, I have the coordinates of both enhancers summits (which could be extended or not) and the coordinates of genes TSSs
    invs_from <- makeGRangesFromDataFrame(invs, seqnames.field = "chr_from", start.field = "chr_from_bkpt", end.field = "chr_from_bkpt", keep.extra.columns = T)
    invs_to <- makeGRangesFromDataFrame(invs, seqnames.field = "chr_to", start.field = "chr_to_bkpt", end.field = "chr_to_bkpt", keep.extra.columns = T)
    
    # Overlap separately from and to, keep which is the loop_id, then merge the 2 datasets keeping this info for both SVs coords 
    hits_from <- findOverlaps(query=invs_from, subject=loops_enh_deg_gr)
    overlaps_from <- cbind(data.frame(mcols(invs_from[queryHits(hits_from)])), data.frame("unique_loop_id_from" =mcols(loops_enh_deg_gr[subjectHits(hits_from)])[,"unique_loop_id"]))
    hits_to <- findOverlaps(invs_to, loops_enh_deg_gr)
    overlaps_to <- cbind(data.frame(mcols(invs_to[queryHits(hits_to)])), data.frame("unique_loop_id_to" =mcols(loops_enh_deg_gr[subjectHits(hits_to)])[,"unique_loop_id"]))
    
    overlaps_join <- full_join(overlaps_from, overlaps_to, by = "sv_id", relationship = "many-to-many")
    
    overlaps_join$ovrlp_from.x[is.na(overlaps_join$ovrlp_from.x)] <- overlaps_join$ovrlp_from.y[is.na(overlaps_join$ovrlp_from.x)]
    overlaps_join$ovrlp_to.x[is.na(overlaps_join$ovrlp_to.x)] <- overlaps_join$ovrlp_to.y[is.na(overlaps_join$ovrlp_to.x)]
    overlaps_join <- overlaps_join %>% dplyr::select(., -c(ovrlp_from.y, ovrlp_to.y)) %>% 
      rename(., c("ovrlp_from" = ovrlp_from.x, "ovrlp_to" = ovrlp_to.x)) %>% 
      relocate(., c(ovrlp_from, ovrlp_to), .after = unique_loop_id_to)
    # check if true - shuold be 0
    sum(!is.na(overlaps_join[overlaps_join$ovrlp_from==F, ]$unique_loop_id_from))
    sum(!is.na(overlaps_join[overlaps_join$ovrlp_to==F, ]$unique_loop_id_to))
    
    #' @Option 1: none of from/to SV coords overlap with a loop
    #' These SVs are not present in overlaps_join, which contains only SVs-loops matches
    
    #' @Option 2: one of the two SV coords (either from or to) overlaps with a loop
    sum(!invs$chr_from_bkpt <= invs$chr_to_bkpt) # always T
    opt_2 <- overlaps_join %>% dplyr::filter(., ovrlp_from+ovrlp_to == 1) 
    # Compute distances: inner bkpt_from - loop_bin2 & loop_bin2 - outer bkpt_to
 
    # inner_dist = bin2(point-in-the-middle) - bkpt1
    opt_2[opt_2$ovrlp_from, ] %>% left_join(., loops_enh_deg, by = c("unique_loop_id_from" = "unique_loop_id")) %>% 
      mutate(., "dist_inner" = ((start2+1000) - chr_from_bkpt) ) %>% View()
    # WARNING, something wrong, alla NAs, missing chr_from_bkpt
  
    # Ideally, we should compute the distance among bkpts and enhancers-summit and DEGs-TSSs
    # But this is complicated, because EACH loop can overlap BOTH 1 enhancer and 1 promoter 
    # Should use unbudled loops to handle this. For now, I will use the loop-bin as a reference 
    
    cat("\n")
  }
    
}


##


# TODO: implement?
# Defining E - P regions - by considering where enhancer falls within bin 
#loops$enh1 <- str_split(loops$name1, ":", simplify = T)[,2] 
#loops <- loops%>% relocate(., enh1, .before = name1)
#enh1_pos <- as.numeric(loops$enh1[1:5]) - loops$start1[1:5]
#if(enh1_pos >= 0 | enh1_pos <= loops_kb*1000) --> enh within bin 1
#if(enh1_pos < 0) --> enh within bin -1
#if(enh1_pos > loops_kb*1000) --> enh within bin +1

#' TODO: 
#' - What about "intrachromosomal rearrangement with non-inverted orientation" ? 
#' - Now we are considering loops connecting enhancers to DEGs promoters. Expand to enhancer - gene?


