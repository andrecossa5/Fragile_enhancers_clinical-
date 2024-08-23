
suppressMessages({

  
library(tidyverse)  # Load tidyverse for data manipulation and visualization

# Set random seed for reproducibility
SEED <- 4321
set.seed(SEED)

# Define parameters
WIN <- 1000  # Window size for considering SNV and enhancer overlap
MARKERS <- c("CtIP", "GRHL")  # List of markers to process
motif_thresh <- 70  # Threshold score percentage for motif matching

# Define file paths for input data
path_variants_anno <- fs::path("/hpcnfs/scratch/P_PGP_FRAGILE_ENHANCERS/results/NEW/data/damaging_variants_annotation/")
path_anno_ehnacers <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/results/integrated/annotated_enhancers/NEW/2kb_GRHL2_enhancers.from_SCR_specific_loops.linked_to_DOWN_DEGs.tsv")


##


# Loop through each marker to process and analyze variant data
for(marker in MARKERS){
  print(paste0("--- Computing fraction of damaging vs. non-damaging variants for ", marker, " enhancers ---"))
  cat("\n")
  
  # Load variant annotation data for the current marker
  variants_anno <- read_tsv(
    fs::path(path_variants_anno, paste0("Table_enh_SSMs_", marker, ".all_overlaps.", WIN, "bp_WIN.with_damaging_anno.motif_thresh_", motif_thresh, ".tsv")))  
  variants_anno$damaging <- factor(variants_anno$damaging, levels = c(1,2,3))
  
  # Print total number of SNVs and the proportion of each damaging level
  print(paste0("Total number of SNVs within all enhancers, across all patients: ", dim(variants_anno)[1]))
  print("How many damaging?")
  print(paste0("Level ", names(table(variants_anno$damaging)), ": ", round(table(variants_anno$damaging) / dim(variants_anno)[1] * 100,2), "%"))
  cat("\n")
  
  # Annotate enhancers with functional status for GRHL marker
  if(marker == "GRHL"){
    # Load enhancer annotation data
    anno_enhancers <- read_tsv(path_anno_ehnacers) %>% mutate(., anno = "anno") 
    variants_anno$enh_anno <- "non-functional"  # Default annotation
    variants_anno[variants_anno$name_enh %in% anno_enhancers$name, ]$enh_anno <- "functional"  # Update to functional where applicable
    variants_anno$enh_anno <- factor(variants_anno$enh_anno, levels = c("functional", "non-functional"))
    
    # Print the number of enhancers in each category
    n_enh <- variants_anno %>% group_by(., enh_anno) %>% dplyr::summarise(., "n" = length(unique(name_enh)))
    print("Number of enhancers x group: ")
    print(paste0(n_enh$enh_anno, ": ", n_enh$n))
    cat("\n")
    
    # Fraction of enhancers with at least 1 damaging mutation in at least 1 patient
    n_enh_1 <- variants_anno %>% dplyr::filter(., damaging == 1) %>% 
      group_by(., enh_anno) %>% dplyr::summarise(., "n" = length(unique(name_enh)))
    n_enh_2 <- variants_anno %>% dplyr::filter(., damaging == 2) %>% 
      group_by(., enh_anno) %>% dplyr::summarise(., "n" = length(unique(name_enh)))
    n_vars <- variants_anno %>% group_by(., enh_anno) %>% dplyr::summarise(., n = n())
    
    # Check that both enhancer categories (functional and non-functional) are present
    if(dim(n_enh_1)[1] == length(unique(variants_anno$enh_anno))){
      print("Fraction of enhancers with at least one damaging mutation - Level 1: ")
      print(paste0(n_enh$enh_anno, ": ", round(n_enh_1$n / n_enh$n,3) * 100, "%"))
      
      # Fraction of damaging mutations over total n. of mutations x enhancer group
      n_vars_1 <- variants_anno %>% dplyr::filter(., damaging == 1) %>% 
        group_by(., enh_anno) %>% dplyr::summarise(., n = n())
      print("Fraction of damaging mutations over tot. mutations - Level 1: ")
      print(paste0(n_vars$enh_anno, ': ', round(n_vars_1$n / n_vars$n,4) * 100, "%"))
    } else {
      print("One level missing. Only one class of enhancers has damaging mutation - Level 1:")
      print(n_enh_1$enh_anno)
    }
    cat("\n")
    
    # Check for presence of damaging mutations (Level 2)
    if(dim(n_enh_2)[1] == length(unique(variants_anno$enh_anno))){
      # Print fraction of enhancers with at least one damaging mutation (Level 2)
      print("Fraction of enhancers with at least one damagin mutation - Level 2: ")
      print(paste0(n_enh$enh_anno, ": ", round(n_enh_2$n / n_enh$n,3) * 100, "%"))  
      
      # Fraction of damaging mutations over total n. of mutations x enhancer group
      n_vars_2 <- variants_anno %>% dplyr::filter(., damaging == 2) %>% 
        group_by(., enh_anno) %>% dplyr::summarise(., n = n())
      print("Fraction of damaging mutations over tot. mutations - Level 2: ")
      print(paste0(n_vars$enh_anno, ': ', round(n_vars_2$n / n_vars$n,4) * 100, "%"))  
    } else {
      print("One level missing. Only one class of enhancers has damaging mutation - Level 2:")
      print(n_enh_2$enh_anno)
    }
    cat("\n")
    
    # Perform statistical tests to compare proportions of damaging mutations in functional vs non-functional enhancers
    print(table(variants_anno$damaging, variants_anno$enh_anno))
    print(fisher.test(table(variants_anno$damaging, variants_anno$enh_anno)))
    
    # Analysis with binary classification (combining Levels 1 and 3)
    variants_anno_1_3 <- variants_anno %>% dplyr::filter(., damaging == 1 | damaging == 3)
    variants_anno_1_3$damaging <- droplevels(variants_anno_1_3$damaging)
    print(table(variants_anno_1_3$damaging, variants_anno_1_3$enh_anno))
    print(fisher.test(table(variants_anno_1_3$damaging, variants_anno_1_3$enh_anno)))
    
    # Binary classification:
    print("Trying with a binary classification")
    # Level 2 becomes damaging (Level 1)
    variants_anno_1_3 <- variants_anno %>% mutate(., damaging = ifelse(damaging == 2, 1, damaging))
    print(table(variants_anno_1_3$damaging, variants_anno_1_3$enh_anno))
    print(fisher.test(table(variants_anno_1_3$damaging, variants_anno_1_3$enh_anno)))
    # Level 2 becomes damaging (Level 1)
    variants_anno_1_3 <- variants_anno %>% mutate(., damaging = ifelse(damaging == 2, 3, damaging))
    print(table(variants_anno_1_3$damaging, variants_anno_1_3$enh_anno))
    print(fisher.test(table(variants_anno_1_3$damaging, variants_anno_1_3$enh_anno)))
    
  }  
}


##

})

##



