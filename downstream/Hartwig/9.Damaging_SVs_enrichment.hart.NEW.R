
suppressMessages({  # Suppresses messages that might clutter the output
  
library(tidyverse)  # Load the tidyverse library for data manipulation and visualization

# Set a seed for reproducibility
SEED <- 4321
set.seed(SEED)

# Define the marker and loop distance variables
marker <- "GRHL"  # Loops annotated with respect to GRHL2 and not CtIP
loops_kb <- 2  # Loops annotated at a *kb resolution

# Define file paths for variant annotations and enhancer data
path_variants_anno <- fs::path(paste0("/hpcnfs/scratch/P_PGP_FRAGILE_ENHANCERS/results/NEW/data/damaging_variants_annotation/SVs/Table_SVs_loops_GRHL.all_overlaps.", loops_kb, "kb_res.with_damaging_anno.tsv"))  
path_anno_ehnacers <- fs::path(paste0("/hpcnfs/scratch/PGP/Ciacci_et_al/results/integrated/NEW/", loops_kb, "kb/data/anno_enhancers/", loops_kb, "kb_GRHL2_enhancers.from_SCR_specific_loops.linked_to_DOWN_DEGs.tsv"))
path_func_loops <- fs::path(paste0("/hpcnfs/scratch/PGP/Ciacci_et_al/results/integrated/NEW/", loops_kb, "kb/data/anno_enhancers/", loops_kb, "kb_functional_loops.aka_SCR_specific_loops_linked_to_DOWN_DEGs.tsv"))


##


# Load the variants annotation data
variants_anno_full <- read_tsv(path_variants_anno)

# Get the unique variant types
types <- c("all", unique(variants_anno_full$variant_type_simple))

# Loop through each variant type
for(type in types){
  cat("\n")  # Print a new line
  print(paste0("--- Computing fraction of damaging vs. non-damaging SVs ( ", type, " ) for ", marker, " enhancers ---"))
  cat("\n")
  
  # Filter the dataset based on the current variant type
  if(type != "all"){
    variants_anno <- variants_anno_full %>% dplyr::filter(., variant_type_simple == type) 
  } else {
    variants_anno <- variants_anno_full
  }
  
  # Print the total number of structural variants (SVs) across all patients
  print(paste0("Total number of SVs (", type, " ) within all E - P regions, across all patients: ", length(unique(variants_anno$sv_id)) ))
  
  # Compute and print the percentage of damaging SVs
  print("How many damaging?")
  dam_vec <- variants_anno %>% group_by(., sv_id) %>% dplyr::summarise(., "how_many_dam" = sum(damaging))  %>% .[,"how_many_dam"]
  print(paste0(round( sum(dam_vec >= 1) / length(unique(variants_anno$sv_id)) * 100), " %"))
  cat("\n")
  
  
  ##
  
  
  # Load functional loops data
  functional_loops <- read_tsv(path_func_loops)
  
  # Annotate loops associated with SVs as functional or non-functional
  variants_anno$loop_anno_from <- "non-functional"
  variants_anno$loop_anno_to <- "non-functional"
  variants_anno[variants_anno$unique_loop_id_from %in% functional_loops$unique_loop_id, ]$loop_anno_from <- "functional"
  variants_anno[variants_anno$unique_loop_id_to %in% functional_loops$unique_loop_id, ]$loop_anno_to <- "functional"
  
  # Create a data frame with unique loop IDs and their annotations
  loops_only <- data.frame(
    "unique_loop_id" = c(variants_anno$unique_loop_id_from, variants_anno$unique_loop_id_to), 
    "loop_anno" = c(variants_anno$loop_anno_from, variants_anno$loop_anno_to), 
    "dam" = c(variants_anno$dam_from, variants_anno$dam_to)
  )
  
  # Remove duplicate and missing entries
  loops_only <- loops_only[!duplicated(loops_only), ] %>% na.omit(.)
  
  # Group by unique loop ID and loop annotation, and summarize any damaging SVs
  loops_only <- loops_only %>% group_by(unique_loop_id, loop_anno) %>%
    dplyr::summarise(any_dam = any(dam), .groups = 'drop')

  # Print the total number of loops overlapping with an SV and the number per category
  print(paste0("Total number of loops overlapping with a SV: ", length(unique(loops_only$unique_loop_id))))
  print("Number of loops overlapping an SV, per category:")
  print(table(loops_only$loop_anno))
  
  
  # Print the number and percentage of loops with at least one damaging SV per category
  print("Number of loops with at least one damaging SV - per group:")
  print(table(loops_only[loops_only$any_dam, ]$loop_anno))
  print("Percentage of loops with at least one damaging SV - per group:")
  print(paste0(round(table(loops_only[loops_only$any_dam, ]$loop_anno) / table(loops_only$loop_anno) * 100,2), " %"))
  
    
  ##
  
  
  # Create a data frame with additional SV information
  loops_only <- data.frame(
    "unique_loop_id" = c(variants_anno$unique_loop_id_from, variants_anno$unique_loop_id_to), 
    "loop_anno" = c(variants_anno$loop_anno_from, variants_anno$loop_anno_to), 
    "dam" = c(variants_anno$dam_from, variants_anno$dam_to), 
    "SV" = c(variants_anno$sv_id, variants_anno$sv_id)
  )
  
  # Remove duplicate and missing entries
  loops_only <- loops_only[!duplicated(loops_only), ] %>% na.omit(.)
  
  # Calculate and print the average number of SVs overlapping with each loop
  avg <- loops_only %>% group_by(unique_loop_id) %>%
    dplyr::summarise(., n_SVs = length(unique(SV)), .groups = "drop") %>%
    .[["n_SVs"]] %>% mean(.) %>% round(.,2)
  print(paste0("Each loop overlaps with an average of: ", avg, " SVs"))
  
  # Calculate and print the percentage of damaging SVs over total SVs per loop category
  print("Percentage of damaging SVs over tot. overlapping SVs - per group: ")
  
  perc_stat <- loops_only %>% group_by(., loop_anno) %>%
    dplyr::summarise(., tot_SVs = length(SV), n_dam = sum(dam)) %>%
    mutate(., perc_dam = n_dam / tot_SVs * 100)
  print(paste0(perc_stat$loop_anno, ": ", round(perc_stat$perc_dam,2), " %"))

}


##

})

##




