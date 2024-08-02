
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

WIN = 1
UNIT = "kb"
MARKER = "CtIP" # or GRHL2

ANALYSIS = "stratified"
GROUP = "BRCA_D"

IN_FOLDER = "/Users/ieo6983/Desktop/fragile_enhancer_clinical/results/ICGC/enhancers_SSMs_overlaps/data/"
if ANALYSIS == "stratified":
    OUT_FOLDER = "/Users/ieo6983/Desktop/fragile_enhancer_clinical/results/ICGC/{}.mutational_signature_analysis.{}{}_WIN.{}/".format(MARKER, WIN, UNIT, GROUP)    
else:
    OUT_FOLDER = "/Users/ieo6983/Desktop/fragile_enhancer_clinical/results/ICGC/{}.mutational_signature_analysis.{}{}_WIN/".format(MARKER, WIN, UNIT)
if not os.path.exists(OUT_FOLDER):
    os.makedirs(OUT_FOLDER)


##


# Read variants 
#SSMs = pd.read_csv(IN_FOLDER+"Table_enh_SSMs_{}.all_overlaps.{}{}_WIN.tsv".format(MARKER, WIN,UNIT), sep = "\t") 
SSMs = pd.read_csv(IN_FOLDER+"Table_enh_SSMs_{}.all_overlaps.{}{}_WIN.{}.tsv".format(MARKER, WIN,UNIT,GROUP), sep = "\t") 

# Define 'high' and 'low' clusters
clusters = {
    "GRHL2" : {
        "high" : ["GRHL_cluster_1_Enh", "GRHL_cluster_2_Enh", "GRHL_cluster_4_Enh"],
        "low" : ["GRHL_cluster_5.2_Enh", "GRHL_cluster_5.3_Enh", "GRHL_cluster_5.1_Enh", "GRHL_cluster_5.0"]
    },
    "CtIP" : {
        "high" : ["CtIP_cluster_2_Enh", "CtIP_cluster_3_Enh", "CtIP_cluster_5_Enh"],
        "low" : ["CtIP_cluster_6.2_Enh", "CtIP_cluster_6.3_Enh", "CtIP_cluster_6.1_Enh", "CtIP_cluster_6.0"]
    }
}


##


# Generate mutational matrices 
for label in clusters[MARKER].keys():
    print(label)

    SSMs_sub = SSMs[SSMs['cluster'].isin(clusters[MARKER][label]) & SSMs['width_sbj'] == 1]

    # Convert to appropriate format 
    format_ssms = {
        'Project' : SSMs_sub['project_code'],
        'Sample' : SSMs_sub['icgc_sample_id'],
        'ID' : SSMs_sub['icgc_donor_id'],
        'Genome' : "GRCh37",
        'mut_type' : "SNP", 
        'chrom' : SSMs_sub['seqnames_sbj'], 
        'pos_start' : SSMs_sub['start_sbj'], 
        'pos_end' : SSMs_sub['end_sbj'], 
        'ref' : SSMs_sub['reference_genome_allele'], 
        'alt' : SSMs_sub['mutated_to_allele'], 
        'Type' : "SOMATIC"
    }
    format_ssms = pd.DataFrame(format_ssms)

    # Save formatted file - required for matrix generation. The folder with input file will correspond to the output folder
    dir_results = OUT_FOLDER+'SNVs_{}_{}/'.format(MARKER, label)
    if not os.path.exists(dir_results):
        os.makedirs(dir_results)
    format_ssms.to_csv(dir_results+'input.SNVs_{}_{}.txt'.format(MARKER,label), index=None, sep = '\t')

    #

    ## Generate mutational matrix - For each sample
    path_to_file = OUT_FOLDER # Full path of the saved input files in the desired output folder.
    
    print("Generating mutational matrix for {} SNVs - Cluster {}".format(MARKER, label))
    matrices = matGen.SigProfilerMatrixGeneratorFunc("SNVs_{}_{}".format(MARKER, label), "GRCh37", path_to_file+'SNVs_{}_{}/'.format(MARKER, label), 
                                                    exome=False, bed_file=None, chrom_based=False, 
                                                    plot=True, tsb_stat=False, seqInfo=True)
    

##



# Use SigProfilerAssignment's main function for mutational signatures assignment
# COSMIC v3.3 mutational signatures are used by default as the input reference signatures

for label in clusters[MARKER].keys():
    file_pattern = "*SBS96.all"

    dir_matrices = path_to_file + "SNVs_{}_{}/output/SBS/".format(MARKER,label)        
    file_list = glob.glob(dir_matrices+file_pattern)[0]

    dir_results_signatures = path_to_file+"SNVs_{}_{}/SigProfilerAssignment_output".format(MARKER, label) 
    if not os.path.exists(dir_results_signatures):
        os.makedirs(dir_results_signatures)

    Analyze.cosmic_fit(samples = file_list, # path to the input somatic mutations file (if using segmentation file/mutational matrix)
                   output = dir_results_signatures, 
                   input_type = "matrix")


##


# Extract activity values for each signature
top_n = 5
for label in clusters[MARKER].keys():
    dir_results_signatures_activities = OUT_FOLDER+"SNVs_{}_{}/SigProfilerAssignment_output/Assignment_Solution/Activities/".format(MARKER, label) 

    SBS_freqs = pd.read_csv(dir_results_signatures_activities+"Assignment_Solution_Activities.txt", 
                            sep = "\t", index_col='Samples')

    # most frequent, aka present in most samples, although not necessarily with more mutations x signature
    total_freqs = SBS_freqs.sum().sort_values(ascending=False)[0:top_n]

    print(f'\n{label} - {MARKER} enhancers\n')

    print("Most active Signatures - Signatures with most mutations")
    print(
        total_freqs
    )

    print("Most frequent Signatures - Signatures present in most samples")
    print(
        (SBS_freqs != 0).sum().sort_values(ascending=False)[0:top_n]
    )


##


## De-novo signature analysis
from SigProfilerExtractor import sigpro as sig # for de-novo extraction of signatures
from SigProfilerAssignment import Analyzer as Analyze # to assign SSMs to a set of de-novo signatures

# Extract signatures de-novo from mutational matrix 
for label in clusters[MARKER].keys():
    file_pattern = "*SBS96.all"

    dir_matrices = OUT_FOLDER + "SNVs_{}_{}/output/SBS/".format(MARKER, label)        
    input_file = glob.glob(dir_matrices+file_pattern)[0]

    out_folder_sig = OUT_FOLDER+"SNVs_{}_{}/de_novo/sig_extraction".format(MARKER, label) 
    if not os.path.exists(out_folder_sig):
        os.makedirs(out_folder_sig)

    sig.sigProfilerExtractor(input_data=input_file, input_type="matrix", output=out_folder_sig, 
                             reference_genome="GRCh37", minimum_signatures=1, maximum_signatures=25)


