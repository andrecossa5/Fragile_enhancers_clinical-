
library(tidyverse)
library(rtracklayer)

### Gene annotation - USCS Table browser - ENSEMBL Genes ###

# Load gene annotation data from a TSV file
anno <- read_tsv("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/others/TSSs_elisa/Gene_anno_hg19_ENSEMBL.tsv")

# Transcripts x gene
# Group the annotation data by gene and count the number of transcripts per gene
tx_x_g <- anno %>% group_by(name2) %>% 
  dplyr::summarise(., n = n())

# Display summary statistics and plot histogram of transcripts per gene
summary(tx_x_g)
hist(table(anno$name2))


# Extract transcription start sites (TSSs) for genes on the positive strand
TSS_plus <- anno[anno$strand == "+", c("name2", "chrom", "strand", "txStart")]
colnames(TSS_plus)[colnames(TSS_plus) == "txStart"] <- "tss"

# Extract transcription start sites (TSSs) for genes on the negative strand
TSS_minus <- anno[anno$strand == "-", c("name2", "chrom", "strand", "txEnd")]
colnames(TSS_minus)[colnames(TSS_minus) == "txEnd"] <- "tss"

# Combine TSSs from both strands into one data frame
TSSs <- rbind(TSS_plus, TSS_minus)

# Genome anno contains also non-standard chromosomes
# Print table of chromosome names to identify non-standard chromosomes
table(TSSs$chrom)

# Select only TSSs located on standard chromosomes (1-22, X, Y)
chroms <- c(paste("chr", seq(1,22), sep=""), "chrX", "chrY")
TSSs <- TSSs[TSSs$chrom %in% chroms, ]
# Verify selection
table(TSSs$chrom)

# Rename 'name2' to 'gene_id' for clarity
colnames(TSSs)[1] <- "gene_id"

# Optionally save the TSSs data to a TSV file
# write_tsv(TSSs, "./TSSs_from_USCS_hg19_EMSEMBL.tsv")

##

# Load TSSs data from a previously saved TSV file
tss <- read_tsv("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/others/TSSs_elisa/TSSs_from_USCS_hg19_EMSEMBL.tsv")
# Display the number of unique gene IDs in the TSSs data
length(unique((tss$gene_id)))
       