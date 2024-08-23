
import pandas as pd
import matplotlib.pyplot as plt
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen # Matrix generation for mutational signatures
from SigProfilerAssignment import Analyzer as Analyze # Signature analysis
import sigProfilerPlotting as sigPlt # Plotting for mutational signatures
import glob
import os

# Install reference genome
#from SigProfilerMatrixGenerator import install as genInstall
#genInstall.install('GRCh37', bash=True)

# Parameters
WIN = 50  # Window size for analysis
UNIT = "bp"  # Unit of measurement for window size
MARKERS = ["CtIP", "GRHL"] # List of markers to analyze

# ANNOTATION: Options for different analyses (commented out)
# ANALYSIS = "stratified"
# GROUP = "BRCA_D"

# Loop through each marker for analysis
for MARKER in MARKERS:
    # Define input and output folders
    IN_FOLDER = "/Users/ieo6983/Desktop/fragile_enhancer_clinical/results/Hartwig/data/"
    OUT_FOLDER = "/Users/ieo6983/Desktop/fragile_enhancer_clinical/results/Hartwig/data/mutational_signature_analysis/{}.mutational_signature_analysis.{}{}_WIN/".format(MARKER, WIN, UNIT)    
    
     # Create output directory if it does not exist
    if not os.path.exists(OUT_FOLDER):
        os.makedirs(OUT_FOLDER)


    ##


    # Read the variants from the input file
    SSMs = pd.read_csv(IN_FOLDER+"SNVs_coords.common_enhancers_mutated_across_hart_and_icgc.{}.tsv".format(MARKER), sep = "\t")
    
    # Process the data: extract chromosome, end position, reference and alternative alleles from the ID column
    SSMs['SAMPLE'] = "unique_sample"
    SSMs['chrom'] = SSMs.ID.str.split(':', expand=True)[[0]]
    SSMs['end'] = SSMs.ID.str.split(':', expand=True)[1].str.split('_', expand=True)[0].astype(int)
    SSMs['ref'] = SSMs.ID.str.split(':', expand=True)[1].str.split('_', expand=True)[1].str.split('/', expand=True)[0]
    SSMs['alt'] = SSMs.ID.str.split(':', expand=True)[1].str.split('_', expand=True)[1].str.split('/', expand=True)[1]

    # Print the status and number of SNVs being analyzed
    print("--- Performing mutational signature analysis --- \n")
    print(f'Regions of interest: {MARKER} enhancers mutated in both Hartwig and ICGC')
    print(f'Total number of SNVs analyzed: {SSMs.shape[0]} \n')

    ##


    # Prepare the data for matrix generation
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

    # Generate the mutational matrix for each sample
    path_to_file = OUT_FOLDER # Path to the saved input files

    print("Generating mutational matrix for {} SNVs".format(MARKER))
    matrices = matGen.SigProfilerMatrixGeneratorFunc("SNVs_{}".format(MARKER), "GRCh37", path_to_file+'SNVs_{}/'.format(MARKER), 
                                                    exome=False, bed_file=None, chrom_based=False, 
                                                    plot=True, tsb_stat=False, seqInfo=True)

        

    ##



    # Use SigProfilerAssignment's main function for mutational signatures assignment
    # COSMIC v3.3 mutational signatures are used by default as the input reference signatures


    file_pattern = "*SBS96.all" # Pattern to match output files

    dir_matrices = path_to_file + "SNVs_{}/output/SBS/".format(MARKER)        
    file_list = glob.glob(dir_matrices+file_pattern)[0]  # Get the file matching the pattern

    dir_results_signatures = path_to_file+"SNVs_{}/SigProfilerAssignment_output".format(MARKER) 
    if not os.path.exists(dir_results_signatures):
        os.makedirs(dir_results_signatures)

    # Run the signature assignment
    Analyze.cosmic_fit(samples = file_list, # Path to the input somatic mutations file
                    output = dir_results_signatures, 
                    input_type = "matrix")


    ##



    # Extract and analyze activity values for each signature
    top_n = 5  # Number of top signatures to display
    dir_results_signatures_activities = OUT_FOLDER+"SNVs_{}/SigProfilerAssignment_output/Assignment_Solution/Activities/".format(MARKER) 

    # Read the activities of the assigned signatures
    SBS_freqs = pd.read_csv(dir_results_signatures_activities+"Assignment_Solution_Activities.txt", 
                            sep = "\t", index_col='Samples')

    # Display most active signatures (those with the most mutations)
    total_freqs = SBS_freqs.sum().sort_values(ascending=False)[0:top_n]

    print(f'\n{MARKER} enhancers\n')

    print("Most active Signatures - Signatures with most mutations")
    print(
        total_freqs
    )

    # Display most frequent signatures (present in most samples)
    print("Most frequent Signatures - Signatures present in most samples")
    print(
        (SBS_freqs != 0).sum().sort_values(ascending=False)[0:top_n]
    )


    ##

