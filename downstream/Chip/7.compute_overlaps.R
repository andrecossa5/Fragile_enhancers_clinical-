
library(tidyverse)
library(GenomicRanges)


##


MARKERS <- c("CtIP", "GRHL")
location <- "local" # 'local' or 'hpc'

path_input_peaks <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/macs2/")
path_peaks_union1 <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/downstream/peaks_union/")
path_output_merged <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHLq05/downstream/peaks_union/")
path_old_enhancers <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/others/")

if(location == "hpc"){
  path_input_peaks <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/bwa/mergedLibrary/macs2/narrowPeak/")
  path_peaks_union1 <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/downstream/peaks_union/")
  path_output_merged <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/downstream/peaks_union/")
  path_old_enhancers <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/others/")
}

# metric1 and metric2 could be pValue and qValue
narrowPeak_colnames <- c("chrom", "start", "end", "name", "score", "strand", 
                         "signal_value", "metric1", "metric2", "summit")
narrowPeak_colnames_nostrand <- c("chrom", "start", "end", "name", "score", 
                                  "signal_value", "metric1", "metric2", "summit")


##


## Check overlap among new and old CtIP / GRHL peaks 

old_ctip <- read_tsv(fs::path(path_old_enhancers, "Cluster_CtIP_Enh_All.txt"), 
                     comment = "#", col_names = c("chrom", "start", "end", "name"))
old_grhl <- read_tsv(fs::path(path_old_enhancers, "Cluster_GRHL_Enh_All.txt"), 
                     comment = "#", col_names = c("chrom", "start", "end", "name"))
old <- list("CtIP" = old_ctip, 
            "GRHL" = old_grhl)

new1 <- list("CtIP" = read_tsv(fs::path(path_peaks_union1, "merged_peaks.hg19_CtIP.filtered_K27ac_k4me1_k4me3.hq_signal.new_summit.narrowPeak"), 
                               col_names = narrowPeak_colnames_nostrand),
             "GRHL" = read_tsv(fs::path(path_peaks_union1, "merged_peaks.hg19_GRHL.filtered_K27ac_k4me1_k4me3.hq_signal.new_summit.narrowPeak"), 
                               col_names = narrowPeak_colnames_nostrand))

old_gr_all <- list()
new1_gr_all <- list()

extend_summits <- T
WIN <- 1000
for(marker in MARKERS){
  cat("\n")
  print(paste0("Computing overlap among ", marker, " old and new peaks"))
  
  old_gr <- old[[marker]] %>% makeGRangesFromDataFrame(., keep.extra.columns=T)
  old_gr_all[[marker]] <- old_gr
  new1_gr <- new1[[marker]] %>% makeGRangesFromDataFrame(., keep.extra.columns=T)
  seqlevels(new1_gr) <- paste0("chr", seqlevels(new1_gr))
  new1_gr_all[[marker]] <- new1_gr

  if(extend_summits ==T){
    start(old_gr) <- start(old_gr) - WIN
    end(old_gr) <- end(old_gr) + WIN
    
    old_gr_all[[marker]]<- old_gr
    old_gr_all[[marker]] <- old_gr
    
    cat("\n")
    print(paste0("Enhancers summits extended of ", WIN, " bp. ", 
                 "Seq length: ", unique(width(old_gr))))
    cat("\n")
  }
  
  ovrlps1 <- !is.na(new1_gr %>% findOverlaps(., old_gr, select = "arbitrary"))
  ovrlps1v <- !is.na(old_gr %>% findOverlaps(., new1_gr, select = "arbitrary"))
  
  print(paste0("Percentage of ", marker, " peaks q05 overlapping with old ones - ", 
        sum(ovrlps1), " over ", length(new1_gr), " : ~", (sum(ovrlps1) / length(new1_gr)) %>% round(.,2) *100, "%"))
  print(paste0("Viceversa - old overlaps with new1 - : ", sum(ovrlps1v), " over ", length(old_gr), " : ", (sum(ovrlps1v) / length(old_gr)) %>% round(.,2) *100, "%"))
}


# Checking also overlap among CtIP and GRHL2 enhancers

# old
forw <- !is.na(findOverlaps(old_gr_all$CtIP, old_gr_all$GRHL, select = "arbitrary"))
backw <- !is.na(findOverlaps(old_gr_all$GRHL, old_gr_all$CtIP, select = "arbitrary"))
print(paste0("OLD CtIP enhancers overlapping with GRHL enhancers: ", 
             sum(forw), " over ", length(forw), " - ", (sum(forw) / length(forw)) %>% round(.,2)*100, "%" ))
print(paste0("OLD GRHL enhancers overlapping with CtIP enhancers: ", 
             sum(backw), " over ", length(backw), " - ", (sum(backw) / length(backw)) %>% round(.,2)*100, "%" ))

# new1 - q-value < 0.05
forw <- !is.na(findOverlaps(new1_gr_all$CtIP, new1_gr_all$GRHL, select = "arbitrary"))
backw <- !is.na(findOverlaps(new1_gr_all$GRHL, new1_gr_all$CtIP, select = "arbitrary"))
print(paste0("NEW1 CtIP enhancers overlapping with GRHL2 enhancers: ", 
             sum(forw), " over ", length(forw), " - ", (sum(forw) / length(forw)) %>% round(.,2)*100, "%" ))
print(paste0("NEW1 GRHL enhancers overlapping with CtIP enhancers: ", 
             sum(backw), " over ", length(backw), " - ", (sum(backw) / length(backw)) %>% round(.,2)*100, "%" ))


##

