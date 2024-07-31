#!/bin/sh
#PBS -l select=1
#PBS -m bae				
#PBS -o /hpcnfs/scratch/PGP/ebugani/fragile_enhancers_clinical/logs/
#PBS -e /hpcnfs/scratch/PGP/ebugani/fragile_enhancers_clinical/errs/

#PBS -N merge_peaks
#PBS -M elisa.bugani@ieo.it		
#PBS -q nocg_workq

source ~/.bashrc
mamba activate deepTools 

IN_DIR='/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/bwa/mergedLibrary/macs2/narrowPeak/'
OUT_DIR='/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/downstream/peaks_union/'


# Merge CtIP peaks
factor='hg19_CtIP'

# Merge all peaks 
#files_to_merge=($(find "$IN_DIR" -type f -name "*.filt.narrowPeak" -printf "%p\n" | grep "$factor"))

# Select only peaks from high-quality signal experiments
files_to_merge=(
	${IN_DIR}'MCF10A_S30627_NT_201120_hg19_CtIP_Total_peaks.filt.narrowPeak'
	${IN_DIR}'MCF10A_S32649_NT_210226_hg19_CtIP_peaks.filt.narrowPeak'
	${IN_DIR}'MCF10A_S40470_EMTnegCTRL2_220211_hg19_CtIP_Tot_peaks.filt.narrowPeak'
	${IN_DIR}'MCF10A_S43554_SCR_220801_hg19_CtIP_Tot_peaks.filt.narrowPeak'
)

cat "${files_to_merge[@]}" | bedtools sort -i - | bedtools merge -i - -c 4,5,6,7,8,9,10 -o collapse -delim "|" > ${OUT_DIR}merged_peaks.${factor}.hq_signal.narrowPeak


# Merge GRHL1/2 peaks
factor='hg19_GRHL'

#files_to_merge=($(find "$IN_DIR" -type f -name "*.filt.narrowPeak" -printf "%p\n" | grep "$factor"))
files_to_merge=(
	${IN_DIR}'MCF10A_S34023_NT_210512_hg19_GRHL2_peaks.filt.narrowPeak'
	${IN_DIR}'MCF10A_S36604_NT_210818_hg19_GRHL1_peaks.filt.narrowPeak'
)

cat "${files_to_merge[@]}" | bedtools sort -i - | bedtools merge -i - -c 4,5,6,7,8,9,10 -o collapse -delim "|" > ${OUT_DIR}merged_peaks.${factor}.hq_signal.narrowPeak








