##For training models
import argparse
from CNN_model_util import preprocessing, eccDNA_CNN_train, eccDNA_CNN_predict, gen_pos_neg_seq
import os    

def Pre_processing(args):
    #cwd is supposed to be DeepCircle 
    cwd = os.getcwd()
    os.system(f"python3 {cwd}/bin/Pre_processing.py \
        --datatype {args.dataType} \
        --boundary {args.boundary} \
        --gap {args.gap} \
        --genome {args.genome} \
        --eccdna {args.eccdna}")
    
def CNN_train(args):
    fname = args.eccDNA_bed.split("/")[-1]
    data_type = fname.split("_")[0]
    arr = preprocessing(ncpus = args.number_of_cpus, genome_fa = args.genome_fasta,\
    genome_gap = args.genome_gap, genome = args.genome, eccdna_bed_dir = args.eccDNA_bed, out_dir = f"{args.output_prediction}/{data_type}")
    eccDNA_CNN_train(*arr,
                    eccdna_bed_dir = args.eccDNA_bed,
                    test_ratio = args.ratio,
                    epochs = args.epochs,
                    batch_size = args.batch_size,
                    output_bed = args.output_prediction)
def CNN_predict(args):
    arr = preprocessing(ncpus = args.number_of_cpus, genome_fa = args.genome_fasta,\
    genome_gap = args.genome_gap, genome = args.genome, eccdna_bed_dir = args.eccDNA_bed, out_dir = f"{args.output_prediction}/{args.chosen_model}")
    eccDNA_CNN_predict(*arr,
                    eccdna_bed_dir = args.eccDNA_bed,
                    chosen_model = args.chosen_model,
                    output_prob = args.output_prediction,
                    genome_fa = args.genome_fasta)
def DNABERT_train(args):
    #cwd is supposed to be DeepCircle 
    cwd = os.getcwd()
    #Data preprocessing
    data_type = args.eccDNA_bed.split("_")[0]
    os.system(f"{cwd}/DNABERT/eccdna/limit_gen.sh -t {data_type}")
    #Fine-tune DNABERT
    os.system(f"python3 {cwd}/DNABERT/examples/run_finetune.py \
    --model_type dnalongcat \
    --tokenizer_name=dna{args.kmer} \
    --model_name_or_path {args.model_path} \
    --task_name dnaprom \
    --do_train \
    --do_eval \
    --data_dir {cwd}/DNABERT/examples/sample_data/ft/eccdna_{data_type}_limit1000/{args.kmer} \
    --max_seq_length 1024 \
    --per_gpu_eval_batch_size=8 \
    --per_gpu_train_batch_size={args.batch_size} \
    --learning_rate 4e-5 \
    --num_train_epochs {args.epochs} \
    --output_dir {cwd}/DNABERT/examples/ft/eccdna_{data_type}_limit1000/{args.kmer} \
    --evaluate_during_training \
    --logging_steps 3000 \
    --save_steps 7000 \
    --warmup_percent 0.1 \
    --hidden_dropout_prob 0.1 \
    --weight_decay 0.001 \
    --n_process {args.number_of_cpus} \
    --overwrite_output_dir \
    --overwrite_cache")
    
def DNABERT_predict(args):
    os.chdir("./DNABERT")
    #Data preprocessing
    data_type = args.eccDNA_bed.split("_")[0]
    os.system(f"./eccdna/limit_gen.sh -t {data_type}")
    #Predict using pre-trained DNABERT
    os.chdir("./examples")
    os.system(f"python3 ./run_finetune.py \
    --model_type dnalongcat \
    --tokenizer_name=dna{args.kmer} \
    --model_name_or_path ./ft/eccdna_{data_type}_limit1000/{args.kmer} \
    --task_name dnaprom \
    --do_predict \
    --data_dir ./sample_data/ft/eccdna_{data_type}_limit1000/{args.kmer} \
    --max_seq_length 1024 \
    --per_gpu_pred_batch_size=128 \
    --output_dir ./ft/eccdna_{data_type}_limit1000/{args.kmer} \
    --predict_dir ./result/eccdna_{data_type}_limit1000/{args.kmer}/ \
    --n_process 48 \
    --overwrite_cache \
    --train_type {args.chosen_model} \
    --test_type {data_type}")
    os.system(f"python3 ./result_transform.py --data {data_type} --len {args.length}")

def infer_motif(args):
    
    if args.method == "DNABERT":
        in_dir = f"./DNABERT/examples/tsv_result/{args.chosen_model}"
        motif_dir = f"{in_dir}/{args.data_type}_{args.chosen_model}_motifs"
        #Generating positive, negative seqs from prediction results
        pred_label_dir = f"{in_dir}/{args.chosen_model}_{args.data_type}_result.tsv"
        pred_seq_dir = f"{in_dir}/{args.data_type}_seq.tsv"
        pos_dir, neg_dir = gen_pos_neg_seq(pred_label_dir, pred_seq_dir)
        
        command = f"~/meme/bin/streme --p {pos_dir} --n {neg_dir} --dna \
        --nmotifs {args.num_motif} --oc {motif_dir}"
        os.system(command)

    if args.method == "CNN":
        in_dir = f"./output/CNN/{args.chosen_model}"
        motif_dir = f"{in_dir}/{args.data_type}_{args.chosen_model}_motifs"
        pos_dir = f"{in_dir}/{args.data_type}_circleseq_eccdna_filt_uniq_pred_positive.fa"
        neg_dir = f"{in_dir}/{args.data_type}_circleseq_eccdna_filt_uniq_pred_negative.fa"
        command = f"~/meme/bin/streme --p {pos_dir} --n {neg_dir} --dna \
        --nmotifs {args.num_motif} --oc {motif_dir}"
        os.system(command)
    
    print(f"Motifs written to {motif_dir}")
# Create the parser
parser = argparse.ArgumentParser(prog = "DeepCircle",
                                 description='Wrapper of CNN and DNABERT within DeepCircle')

subparsers = parser.add_subparsers()

parser_pre_processing = subparsers.add_parser('Pre-processing')
parser_pre_processing.set_defaults(func = Pre_processing)
parser_CNN_train = subparsers.add_parser('CNN-train')
parser_CNN_train.set_defaults(func = CNN_train)
parser_CNN_predict = subparsers.add_parser('CNN-predict')
parser_CNN_predict.set_defaults(func = CNN_predict)
parser_DNABERT_train = subparsers.add_parser('DNABERT-train')
parser_DNABERT_train.set_defaults(func = DNABERT_train)
parser_DNABERT_predict = subparsers.add_parser('DNABERT-predict')
parser_DNABERT_predict.set_defaults(func = DNABERT_predict)
parser_interpret = subparsers.add_parser('infer-motif')
parser_interpret.set_defaults(func = infer_motif)

# Add the arguments

#Parser for pre-processing
parser_pre_processing.add_argument("-ext",
                    "--extend",
                    type=int,
                    default=512,
                    help="Sequence length you want to extend from center"
)

parser_pre_processing.add_argument("-win",
                    "--window",
                    type=int,
                    default=1000,
                    help="Window size of sequence in CNN"
)

parser_pre_processing.add_argument("-l",
                    "--limit",
                    type=int,
                    default=1000,
                    help="Limit of sequence length"
)

parser_pre_processing.add_argument("-type",
                    "--dataType",
                    type=str,
                    help="Species name of eccdna"
)

parser_pre_processing.add_argument("-bound",
                    "--boundary",
                    type=str,
                    help="Genome boundary file name"
)

parser_pre_processing.add_argument("-g",
                    "--gap",
                    type=str,
                    help="Bedfile name of genome gap"
)

parser_pre_processing.add_argument("-gen",
                    "--genome",
                    type=str,
                    help="Genome reference file name"
)

parser_pre_processing.add_argument("-ecc",
                    "--eccdna",
                    type=str,
                    help="Eccdna bedfile name"
)

#Parser for training CNN

parser_CNN_train.add_argument("-c",
                    "--number-of-cpus",
                    type=int,
                    default = 8,
                    help="Number of cpus [default: 8]")

parser_CNN_train.add_argument("-fa",
                    "--genome-fasta",
                    type=str,
                    default = "./DNABERT/eccdna/genome/human/hg38_noalt.fa",
                    help='Directory to the fasta file of hg38 genome\
                    [default: ./DNABERT/eccdna/genome/human/hg38_noalt.fa]')

parser_CNN_train.add_argument("-gap",
                    "--genome-gap",
                    type=str,
                    default = "./DNABERT/eccdna/genome/human/hg38_noalt_gap.bed",
                    help='Directory to the bed file of gaps of hg38 genome\
                    [default: ./DNABERT/eccdna/genome/human/hg38_noalt_gap.bed]')

parser_CNN_train.add_argument('-g',
                    '--genome',
                    type=str,
                    default = "./DNABERT/eccdna/genome/human/hg38_noalt.fa.genome",
                    help="Directory to the file containing lengths of each chromosome in hg38 genome\
                    [default: ./DNABERT/eccdna/genome/human/hg38_noalt.fa.genome]")

parser_CNN_train.add_argument('-bed',
                    '--eccDNA-bed',
                    type=str,
                    default = "./DNABERT/eccdna/db/human/PC-3_circleseq_eccdna_filt_uniq.bed",
                    required = True,
                    help='Directory to the bed file containing eccDNAs\
                    [default: ./DNABERT/eccdna/db/human/PC-3_circleseq_eccdna_filt_uniq.bed]')

parser_CNN_train.add_argument('-r',
                    '--ratio',
                    default = 0.2,
                    type=float,
                    help='Ratio of data used for testing [default: 0.2]')

parser_CNN_train.add_argument('-epo',
                    '--epochs',
                    default = 5,
                    type=int,
                    help='Number of epochs for training [default: 5]')

parser_CNN_train.add_argument('-batch',
                    '--batch-size',
                    type=int,
                    default = 32,
                    help='Number of batch size for training [default: 32]')

parser_CNN_train.add_argument("-dir",
                    "--output-prediction",
                    type = str,
                    default = "./output/CNN",
                    help="Directory to output prediction results of the model [default: ./output/CNN]")

#Parser for CNN predicting
parser_CNN_predict.add_argument("-c",
                    "--number-of-cpus",
                    type=int,
                    default = 8,
                    help="Number of cpus [default: 8]")

parser_CNN_predict.add_argument("-fa",
                    "--genome-fasta",
                    type=str,
                    default = "./DNABERT/eccdna/genome/human/hg38_noalt.fa",
                    help='Directory to the fasta file of hg38 genome [default: ./DNABERT/eccdna/genome/human/hg38_noalt.fa]')

parser_CNN_predict.add_argument("-gap",
                    "--genome-gap",
                    type=str,
                    default = "./DNABERT/eccdna/genome/human/hg38_noalt_gap.bed",
                    help='Directory to the bed file of gaps of hg38 genome\
                    [default: ./DNABERT/eccdna/genome/human/hg38_noalt_gap.bed]')

parser_CNN_predict.add_argument('-g',
                    '--genome',
                    type=str,
                    default = "./DNABERT/eccdna/genome/human/hg38_noalt.fa.genome",
                    help="Directory to the file containing lengths of each chromosome in hg38 genome\
                    [default: ./DNABERT/eccdna/genome/human/hg38_noalt.fa.genome]")

parser_CNN_predict.add_argument('-bed',
                    '--eccDNA-bed',
                    type=str,
                    default = "./DNABERT/eccdna/db/human/PC-3_circleseq_eccdna_filt_uniq.bed",
                    required = True,
                    help='Directory to the bed file containing eccDNAs\
                    [default: ./DNABERT/eccdna/db/human/PC-3_circleseq_eccdna_filt_uniq.bed]')

parser_CNN_predict.add_argument("-dir",
                    "--output-prediction",
                    type = str,
                    default = "./output/CNN",
                    help="Directory to output prediction results of the model [default: ./output/CNN]")

parser_CNN_predict.add_argument('-model',
                    '--chosen-model',
                    choices = ["C4-2", "ES2", "HeLaS3", "LnCap", "OVCAR8", "PC-3", "U937", "Leukocytes", \
                    "Muscle", "Lung-normal", "Lung-tumor"],
                    required = True,
                    type=str,
                    help='Choose pre-trained models for predicting eccDNAs [default: DeepCircle/pre_trained_models/]')

#Parser for training DNABERT
parser_DNABERT_train.add_argument("-k",
                    "--kmer",
                    type=int,
                    default = 6,
                    help='Specify k-mer used to train DNABERT [default: 6]')

parser_DNABERT_train.add_argument("-model_path",
                    "--model_path",
                    type=str,
                    default = "./DNABERT/6-new-12w-0",
                    help='Directory to the pre-trained DNABERT, \
                    [default: ./DNABERT/6-new-12w-0]')

parser_DNABERT_train.add_argument('-batch',
                    '--batch-size',
                    type=int,
                    default = 8,
                    help='Number of batch size for training')

parser_DNABERT_train.add_argument('-epo',
                    '--epochs',
                    default = 2,
                    type=int,
                    help='Number of epochs for training')

parser_DNABERT_train.add_argument("-c",
                    "--number-of-cpus",
                    type=int,
                    default = 2,
                    help="Number of cpus")

#Parser for predicting using DNABERT
parser_DNABERT_predict.add_argument("-k",
                    "--kmer",
                    type=int,
                    default = 6,
                    help='Specify k-mer used to train DNABERT [default: 6]')

parser_DNABERT_predict.add_argument('-batch',
                    '--batch-size',
                    type=int,
                    default = 8,
                    help='Number of batch size for training')

parser_DNABERT_predict.add_argument('-epo',
                    '--epochs',
                    default = 2,
                    type=int,
                    help='Number of epochs for training')

parser_DNABERT_predict.add_argument("-c",
                    "--number-of-cpus",
                    type=int,
                    default = 2,
                    help="Number of cpus")

parser_DNABERT_predict.add_argument('-model',
                    '--chosen-model',
                    choices = ["C4-2", "ES2", "HeLaS3", "LnCap", "OVCAR8", "PC-3", "U937", "Leukocytes", \
                    "Muscle", "Lung-normal", "Lung-tumor"],
                    type=str,
                    help='Choose pre-trained models for predicting eccDNAs [default: DeepCircle/DNABERT/examples/ft/]')

parser_DNABERT_predict.add_argument('-len',
                    '--length',
                    choices = [512, 1024],
                    type=int,
                    default=1024,
                    help='Sequence length of data set')

#Parser for inferring motifs
parser_interpret.add_argument('-method',
                    choices = ["CNN", "DNABERT"],
                    type=str,
                    default = "CNN",
                    help='Choose type of model to that was used to predict eccDNA [default: CNN]')
parser_interpret.add_argument('-data',
                    '--data-type',
                    choices = ["C4-2", "ES2", "HeLaS3", "LnCap", "OVCAR8", "PC-3", "U937", "Leukocytes", \
                    "Muscle", "Lung-normal", "Lung-tumor"],
                    required=True,
                    type=str,
                    help='Choose type of data that was used to predict eccDNA')
parser_interpret.add_argument('-n',
                    '--num-motif',
                    type=int,
                    default=10,
                    help='Type number of motifs to infer [default: 10]')
parser_interpret.add_argument('-model',
                    '--chosen-model',
                    choices = ["C4-2", "ES2", "HeLaS3", "LnCap", "OVCAR8", "PC-3", "U937", "Leukocytes", \
                    "Muscle", "Lung-normal", "Lung-tumor"],
                    type=str,
                    help='Choose pre-trained models for predicting eccDNAs')


if __name__ == "__main__":
    args = parser.parse_args()
    args.func(args)
        



"""

"""


