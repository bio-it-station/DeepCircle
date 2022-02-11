import argparse
from CNN_model_util import *

# Create the parser
parser = argparse.ArgumentParser(prog = "DeepCircle-predict",
                                 usage = "DeepCircle_predict.py [-c <int>] [-fa <dir>] [-g <dir>] [-bed <dir>] [-out] [-model {C4-2, ES2, LnCap, OVCAR8, PC-3, U937, leukocytes, muscle, pool_LCN, pool_LCT}]",
                                 description='Predict the probability of eccDNAs present in the given genomic intervals.')

# Add the arguments
parser.add_argument('-c',
                    '--number-of-cpus',
                    type=int,
                    help='Number of cpus')
parser.add_argument('-fa',
                    '--genome-fasta',
                    type=str,
                    help='Directory to the fasta file of hg38 genome')
parser.add_argument('-g',
                    '--genome',
                    type=str,
                    help='Directory to the file containing lengths of each chromosome in hg38 genome')
parser.add_argument('-bed',
                    '--eccDNA-bed',
                    type=str,
                    help='Directory to the bed file containing eccDNAs')
parser.add_argument('-model',
                    '--chosen-model',
                    choices = ["C4-2", "ES2", "HeLaS3", "LnCap", "OVCAR8", "PC-3", "U937", "leukocytes", "muscle", "pool_LCN", "pool_LCT"],
                    type=str,
                    help='Choose pre-trained models for predicting eccDNAs.')
parser.add_argument('-out',
                    '--output-prob',
                    action = "store_true",
                    help='Output probability of eccDNAs present in given genomic windoes.')

args = parser.parse_args()
arr = preprocessing(ncpus = args.number_of_cpus, genome_fa = args.genome_fasta, genome = args.genome, eccdna_bed_dir = args.eccDNA_bed)
eccDNA_CNN_predict(*arr,
                   eccdna_bed_dir = args.eccDNA_bed,
                   chosen_model = args.chosen_model,
                   output_prob = args.output_prob)
