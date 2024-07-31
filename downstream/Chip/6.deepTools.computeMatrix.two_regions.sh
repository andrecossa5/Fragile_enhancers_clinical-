#!/bin/sh
#PBS -l select=1
#PBS -m bae				
#PBS -o /hpcnfs/scratch/PGP/ebugani/fragile_enhancers_clinical/logs/
#PBS -e /hpcnfs/scratch/PGP/ebugani/fragile_enhancers_clinical/errs/
#PBS -N computeMatrix
#PBS -M elisa.bugani@ieo.it		
#PBS -q nocg_workq

source ~/.bashrc
mamba activate deepTools

IN_DIR='/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/downstream/deepTools/'

factor1='CtIP'
factor2='GRHL'
factor_regions1='CtIP'
factor_regions2='GRHL'
OP='sum' # 'mean' or 'sum' - based on what was used to merge bigWig files.

REGIONS_DIR='/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/downstream/peaks_union/'
REGIONS_FILE1=$(ls ${REGIONS_DIR}*${factor_regions1}.filtered*hq_signal.new_summit.narrowPeak)
REGIONS_FILE2=$(ls ${REGIONS_DIR}*${factor_regions2}.filtered*hq_signal.new_summit.narrowPeak)

BIGWIG_DIR='/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/downstream/deepTools/'
BIGWIG_FILE1=$(ls ${BIGWIG_DIR}'combined_all.'${factor1}'.hq_signal.'${OP}'.bw')
BIGWIG_FILE2=$(ls ${BIGWIG_DIR}'combined_all.'${factor2}'.hq_signal.'${OP}'.bw')
 
cd $IN_DIR

computeMatrix reference-point \
	-S ${BIGWIG_FILE1} ${BIGWIG_FILE2} \
	-R ${REGIONS_FILE1} ${REGIONS_FILE2} \
	--referencePoint center \
	--beforeRegionStartLength 3000 \
	--afterRegionStartLength 3000 \
	--skipZeros \
	-o matrix.CtIP_and_GRHL_bw.CtIP_and_GRHL_regions.hq_signal.from_mean.gz \
	--outFileSortedRegions matrix.CtIP_and_GRHL_bw.CtIP_and_GRHL_regions.hq_signal.from_mean.sortedRegions.bed






