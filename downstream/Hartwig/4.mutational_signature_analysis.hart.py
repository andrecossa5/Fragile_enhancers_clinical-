
import pandas as pd
import matplotlib.pyplot as plt
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen # matrix generation
from SigProfilerAssignment import Analyzer as Analyze # signature analysis
import sigProfilerPlotting as sigPlt # signatures plots
import glob
import os

# Install reference genome
#from SigProfilerMatrixGenerator import install as genInstall
#genInstall.install('GRCh37', bash=True)

WIN = 50
UNIT = "bp"
MARKERS = ["CtIP", "GRHL"] # or GRHL2

#ANALYSIS = "stratified"
#GROUP = "BRCA_D"

for MARKER in MARKERS:

    IN_FOLDER = "/Users/ieo6983/Desktop/fragile_enhancer_clinical/results/Hartwig/data/"
    OUT_FOLDER = "/Users/ieo6983/Desktop/fragile_enhancer_clinical/results/Hartwig/data/mutational_signature_analysis/{}.mutational_signature_analysis.{}{}_WIN/".format(MARKER, WIN, UNIT)    
    if not os.path.exists(OUT_FOLDER):
        os.makedirs(OUT_FOLDER)


    ##


    # Read variants 
    SSMs = pd.read_csv(IN_FOLDER+"SNVs_coords.common_enhancers_mutated_across_hart_and_icgc.{}.tsv".format(MARKER), sep = "\t")
    SSMs['SAMPLE'] = "unique_sample"
    SSMs['chrom'] = SSMs.ID.str.split(':', expand=True)[[0]]
    SSMs['end'] = SSMs.ID.str.split(':', expand=True)[1].str.split('_', expand=True)[0].astype(int)
    SSMs['ref'] = SSMs.ID.str.split(':', expand=True)[1].str.split('_', expand=True)[1].str.split('/', expand=True)[0]
    SSMs['alt'] = SSMs.ID.str.split(':', expand=True)[1].str.split('_', expand=True)[1].str.split('/', expand=True)[1]

    # TODO: add print() of how many SNVs are being analyzed and what they are 
    print("--- Performing mutational signature analysis --- \n")
    print(f'Regions of interest: {MARKER} enhancers mutated in both Hartwig and ICGC')
    print(f'Total number of SNVs analyzed: {SSMs.shape[0]} \n')


    ##


    # Generate mutational matrices 
    format_ssms = {
        'Project' : ".",
        'Sample' : SSMs['SAMPLE'],
        'ID' : ".",
        'Genome' : "GRCh37",
        'mut_type' : "SNP", 
        'chrom' : SSMs['chrom'], 
        'pos_start' : SSMs['end'], 
        'pos_end' : SSMs['end'], 
        'ref' : SSMs['ref'], 
        'alt' : SSMs['alt'], 
        'Type' : "SOMATIC"
    }
    format_ssms = pd.DataFrame(format_ssms)

    # Save formatted file - required for matrix generation. The folder with input file will correspond to the output folder
    dir_results = OUT_FOLDER+'SNVs_{}/'.format(MARKER)
    if not os.path.exists(dir_results):
        os.makedirs(dir_results)
    format_ssms.to_csv(dir_results+'input.SNVs_{}.txt'.format(MARKER), index=None, sep = '\t')

    #

    ## Generate mutational matrix - For each sample
    path_to_file = OUT_FOLDER # Full path of the saved input files in the desired output folder.

    print("Generating mutational matrix for {} SNVs".format(MARKER))
    matrices = matGen.SigProfilerMatrixGeneratorFunc("SNVs_{}".format(MARKER), "GRCh37", path_to_file+'SNVs_{}/'.format(MARKER), 
                                                    exome=False, bed_file=None, chrom_based=False, 
                                                    plot=True, tsb_stat=False, seqInfo=True)
        

    ##



    # Use SigProfilerAssignment's main function for mutational signatures assignment
    # COSMIC v3.3 mutational signatures are used by default as the input reference signatures


    file_pattern = "*SBS96.all"

    dir_matrices = path_to_file + "SNVs_{}/output/SBS/".format(MARKER)        
    file_list = glob.glob(dir_matrices+file_pattern)[0]

    dir_results_signatures = path_to_file+"SNVs_{}/SigProfilerAssignment_output".format(MARKER) 
    if not os.path.exists(dir_results_signatures):
        os.makedirs(dir_results_signatures)

    Analyze.cosmic_fit(samples = file_list, # path to the input somatic mutations file (if using segmentation file/mutational matrix)
                    output = dir_results_signatures, 
                    input_type = "matrix")


    ##


    # Extract activity values for each signature
    top_n = 5
    dir_results_signatures_activities = OUT_FOLDER+"SNVs_{}/SigProfilerAssignment_output/Assignment_Solution/Activities/".format(MARKER) 

    SBS_freqs = pd.read_csv(dir_results_signatures_activities+"Assignment_Solution_Activities.txt", 
                            sep = "\t", index_col='Samples')

    # most frequent, aka present in most samples, although not necessarily with more mutations x signature
    total_freqs = SBS_freqs.sum().sort_values(ascending=False)[0:top_n]

    print(f'\n{MARKER} enhancers\n')

    print("Most active Signatures - Signatures with most mutations")
    print(
        total_freqs
    )

    print("Most frequent Signatures - Signatures present in most samples")
    print(
        (SBS_freqs != 0).sum().sort_values(ascending=False)[0:top_n]
    )


    ##

