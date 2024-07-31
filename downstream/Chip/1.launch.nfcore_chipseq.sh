#!/bin/sh
#PBS -l select=1
#PBS -m bae				
#PBS -o /hpcnfs/scratch/PGP/ebugani/fragile_enhancers_clinical/logs/
#PBS -e /hpcnfs/scratch/PGP/ebugani/fragile_enhancers_clinical/errs/

#PBS -N nf_chip_seq
#PBS -M elisa.bugani@ieo.it		
#PBS -q nocg_workq

source ~/.bashrc
mamba activate nf_env

IN_DIR='/hpcnfs/scratch/PGP/ebugani/fragile_enhancers_clinical/run_nf/'
OUT_DIR='/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/Chip/Chip_for_clusters/results/CtIP_GRHL_q05/'
INPUT_SHEET='/hpcnfs/scratch/PGP/ebugani/fragile_enhancers_clinical/chip_marta_samplesheet.csv'

if [ ! -d $OUT_DIR ]; then
	mkdir -p $OUT_DIR
fi 

cd $IN_DIR

nextflow run nf-core/chipseq -r 2.0.0 \
    --input $INPUT_SHEET \
    --outdir $OUT_DIR \
    --genome GRCh37 \  
    --macs_gsize 2864785220 \
    --narrow_peak \
    -profile singularity \
    -with-tower \
    --macs_pvalue 0.01 
   
# blacklisted regions for the specified genome are bundled with the pipeline 








