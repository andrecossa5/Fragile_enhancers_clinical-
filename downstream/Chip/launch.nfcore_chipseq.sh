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
OUT_DIR='/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/Chip/Chip_for_clusters/results/'
INPUT_SHEET='/hpcnfs/scratch/PGP/ebugani/fragile_enhancers_clinical/chip_marta_samplesheet.csv'

cd $IN_DIR

#export TOWER_ACCESS_TOKEN=eyJ0aWQiOiA5NjQ2fS44NjhiOTUzN2FhY2YxNGFlNDVmOTA3YmY0OWIxYzA4YTJiMmI1ZmZh

nextflow run nf-core/chipseq -r 2.0.0 \
    --input $INPUT_SHEET \
    --outdir $OUT_DIR \
    --genome GRCh37 \  
    --macs_gsize 2864785220 \
    --narrow_peak \
    -profile singularity \
    -with-tower \
    -resume

# blacklisted regions for the specified genome are bundled with the pipeline 








