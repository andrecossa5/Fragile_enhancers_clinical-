#!/bin/sh
#PBS -l select=1
#PBS -m bae                             
#PBS -o /hpcnfs/scratch/PGP/ebugani/fragile_enhancers_clinical/logs/
#PBS -e /hpcnfs/scratch/PGP/ebugani/fragile_enhancers_clinical/errs/
#PBS -N bigwigCompare
#PBS -M elisa.bugani@ieo.it             
#PBS -q nocg_workq

source ~/.bashrc
mamba activate deepTools

IN_DIR='/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/downstream/deepTools/'
OUT_DIR='/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/downstream/deepTools/'

BIGWIG_DIR='/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/bwa/mergedLibrary/bigwig/'
#BIGWIG_FILES=$(ls ${BIGWIG_DIR}*${factor}*)

cd $BIGWIG_DIR

factor='CtIP'
bigwigCompare -b1 MCF10A_S30627_NT_201120_hg19_CtIP_Total.bigWig -b2 MCF10A_S32649_NT_210226_hg19_CtIP.bigWig --operation mean -o ${OUT_DIR}combined_1_2.bw
bigwigCompare -b1 ${OUT_DIR}combined_1_2.bw -b2 MCF10A_S40470_EMTnegCTRL2_220211_hg19_CtIP_Tot.bigWig --operation mean -o ${OUT_DIR}combined_1_to_3.bw
bigwigCompare -b1 ${OUT_DIR}combined_1_to_3.bw -b2 MCF10A_S43554_SCR_220801_hg19_CtIP_Tot.bigWig --operation mean -o ${OUT_DIR}combined_all.${factor}.hq_signal.mean.bw

bigwigCompare -b1 MCF10A_S30627_NT_201120_hg19_CtIP_Total.bigWig -b2 MCF10A_S32649_NT_210226_hg19_CtIP.bigWig --operation add -o ${OUT_DIR}combined_1_2.bw
bigwigCompare -b1 ${OUT_DIR}combined_1_2.bw -b2 MCF10A_S40470_EMTnegCTRL2_220211_hg19_CtIP_Tot.bigWig --operation add -o ${OUT_DIR}combined_1_to_3.bw
bigwigCompare -b1 ${OUT_DIR}combined_1_to_3.bw -b2 MCF10A_S43554_SCR_220801_hg19_CtIP_Tot.bigWig --operation add -o ${OUT_DIR}combined_all.${factor}.hq_signal.sum.bw

factor='GRHL'
bigwigCompare -b1 MCF10A_S34023_NT_210512_hg19_GRHL2.bigWig -b2 MCF10A_S36604_NT_210818_hg19_GRHL1.bigWig --operation mean -o ${OUT_DIR}combined_all.${factor}.hq_signal.mean.bw

bigwigCompare -b1 MCF10A_S34023_NT_210512_hg19_GRHL2.bigWig -b2 MCF10A_S36604_NT_210818_hg19_GRHL1.bigWig --operation add -o ${OUT_DIR}combined_all.${factor}.hq_signal.sum.bw



