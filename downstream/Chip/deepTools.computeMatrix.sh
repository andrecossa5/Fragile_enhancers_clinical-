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

factor='CtIP'
factor_regions='CtIP'
OP='sum' # 'mean' or 'sum' - based on what was used to merge bigWig files.

REGIONS_DIR='/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/downstream/peaks_union/'
REGIONS_FILE=$(ls ${REGIONS_DIR}*${factor_regions}.filtered*hq_signal.new_summit.narrowPeak)

#BIGWIG_DIR='/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/bwa/mergedLibrary/bigwig/'
BIGWIG_DIR='/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/downstream/deepTools/'
BIGWIG_FILE=$(ls ${BIGWIG_DIR}'combined_all.'${factor}'.hq_signal.'${OP}'.bw')

 
cd $IN_DIR

computeMatrix reference-point \
	-R $REGIONS_FILE \	
	-S $BIGWIG_FILE \
	--referencePoint center \
	-o matrix.${factor}_bw.hq_signal.regions_${factor_regions}_peaks_union.from_${OP}.gz \
	--outFileSortedRegions matrix.${factor}_bw.hq_signal.regions_${factor_regions}_peaks_union.from_${OP}.sortedRegions.bed \
	--beforeRegionStartLength 3000 \
	--afterRegionStartLength 3000 \
	--skipZeros 









