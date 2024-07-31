#!/bin/sh
#PBS -l select=1
#PBS -m bae				
#PBS -o /hpcnfs/scratch/PGP/ebugani/fragile_enhancers_clinical/logs/
#PBS -e /hpcnfs/scratch/PGP/ebugani/fragile_enhancers_clinical/errs/
#PBS -N deepTools
#PBS -M elisa.bugani@ieo.it		
#PBS -q nocg_workq

source ~/.bashrc
mamba activate deepTools

IN_DIR='/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL/downstream/deepTools/'

factor='CtIP'
REGIONS_DIR='/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL/downstream/peaks_union/'
REGIONS_FILE=$(ls ${REGIONS_DIR}/*${factor}.filtered*.narrowPeak)
BIGWIG_DIR='/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL/bwa/mergedLibrary/bigwig/'
BIGWIG_FILES=$(ls ${BIGWIG_DIR}/*${factor}*)
 
cd $IN_DIR

#plotHeatmap -m matrix.${factor}.gz \
#	-out heatmap.${factor}.png

plotHeatmap -m matrix.GRHL_bw.MCF10A_S34023.regions_CtIP_peaks_union.gz \
	-out heatmap.GRHL_bw.MCF10A_S34023.regions_CtIP_peaks_union.png


mamba deactivate


