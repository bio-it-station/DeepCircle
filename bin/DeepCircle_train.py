##For training models
import argparse
from CNN_model_util import *

# Create the parser
parser = argparse.ArgumentParser(prog = "DeepCircle-train",
                                 usage = "DeepCircle-train [-c <int>] [-fa <dir>] [-g <dir>] [-bed <dir>] [-r <float>] [-epo <int>] [-batch] [-out]",
                                 description='Train CNN based on provided eccDNAs.')

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
parser.add_argument('-r',
                    '--ratio',
                    type=float,
                    help='Ratio of data used for testing')
parser.add_argument('-epo',
                    '--epochs',
                    type=str,
                    help='Number of epochs for training')
parser.add_argument('-batch',
                    '--batch-size',
                    type=int,
                    help='Number of batch size for training')
parser.add_argument('-dir',
                    '--output-prediction',
                    type=str,
                    help='Directory to output prediction results of the model [default: ./output]')

args = parser.parse_args()
arr = preprocessing(ncpus = args.number_of_cpus, genome_fa = args.genome_fasta, genome = args.genome, eccdna_bed_dir = args.eccDNA_bed)
eccDNA_CNN_train(*arr,
                 eccdna_bed_dir = args.eccDNA_bed,
                 test_ratio = args.ratio,
                 epochs = args.epochs,
                 batch_size = args.batch_size,
                 output_bed = args.output_prediction)


