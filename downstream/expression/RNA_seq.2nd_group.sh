#!/bin/sh
#PBS -l select=1	
#PBS -m ae				
#PBS -o /hpcnfs/scratch/PGP/ebugani/expression/logs/
#PBS -e /hpcnfs/scratch/PGP/ebugani/expression/errs/

#PBS -N RNA_seq_nf.2nd_group
#PBS -M elisa.bugani@ieo.it		
#PBS -q nocg_workq		

source ~/.bashrc
mamba activate nf_env

IN_DIR='/hpcnfs/scratch/PGP/ebugani/expression/data/fasterq_2nd_group/'
OUT_DIR='/hpcnfs/scratch/PGP/ebugani/expression/results_2nd_group/'

cd $IN_DIR

nextflow run nf-core/rnaseq -r 3.14.0 \
    --input /hpcnfs/scratch/PGP/ebugani/expression/samplesheet.no_bom.2nd_group.csv \
    --outdir $OUT_DIR \
    --genome GRCh37 \
    -profile singularity
