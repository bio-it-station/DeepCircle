# DeepCircle
Predict eccDNA based on CNN

# Install

```sh
conda env create -f environment.yml
```
# Quick start
Predict eccDNAs based on given bed files
```Python
python3 DeepCircle_predict.py -c 8 -fa ../genome/hg38_noalt.fa -g ../genome/hg38_noalt.fa.genome -bed ../example/leukocytes_circleseq_eccdna_filt_uniq.bed -model HeLaS3 -out
```

Train CNN based on provided eccDNA
```Python
python3 DeepCircle_train.py -c 8 -fa ../genome/hg38_noalt.fa -g ../genome/hg38_noalt.fa.genome -bed ../example/leukocytes_circleseq_eccdna_filt_uniq.bed -r 0.2 -epo 5 -batch 100 -out
```

# Advance usage
# Training
```
usage: DeepCircle-train [-c <int>] [-fa <dir>] [-g <dir>] [-bed <dir>] [-r <float>] [-epo <int>] [-batch] [-out]

Train CNN based on provided eccDNAs.

optional arguments:
  -h, --help            show this help message and exit
  -c NUMBER_OF_CPUS, --number-of-cpus NUMBER_OF_CPUS
                        Number of cpus
  -fa GENOME_FASTA, --genome-fasta GENOME_FASTA
                        Directory to the fasta file of hg38 genome
  -g GENOME, --genome GENOME
                        Directory to the file containing lengths of each
                        chromosome in hg38 genome
  -bed ECCDNA_BED, --eccDNA-bed ECCDNA_BED
                        Directory to the bed file containing eccDNAs
  -r RATIO, --ratio RATIO
                        Ratio of data used for testing
  -epo EPOCHS, --epochs EPOCHS
                        Number of epochs for training
  -batch BATCH_SIZE, --batch-size BATCH_SIZE
                        Number of batch size for training
  -out, --output-prediction
                        Output prediction results of the model

```

# Prediction
```
usage: DeepCircle_predict.py [-c <int>] [-fa <dir>] [-g <dir>] [-bed <dir>] [-out] [-model {C4-2, ES2, LnCap, OVCAR8, PC-3, U937, leukocytes, muscle, pool_LCN, pool_LCT}]
Predict the probability of eccDNAs present in the given genomic intervals.

optional arguments:
  -h, --help            show this help message and exit
  -c NUMBER_OF_CPUS, --number-of-cpus NUMBER_OF_CPUS
                        Number of cpus
  -fa GENOME_FASTA, --genome-fasta GENOME_FASTA
                        Directory to the fasta file of hg38 genome
  -g GENOME, --genome GENOME
                        Directory to the file containing lengths of each
                        chromosome in hg38 genome
  -bed ECCDNA_BED, --eccDNA-bed ECCDNA_BED
                        Directory to the bed file containing eccDNAs
  -model {C4-2,ES2,HeLaS3,LnCap,OVCAR8,PC-3,U937,leukocytes,muscle,pool_LCN,pool_LCT}, --chosen-model {C4-2,ES2,HeLaS3,LnCap,OVCAR8,PC-3,U937,leukocytes,muscle,pool_LCN,pool_LCT}
                        Choose pre-trained models for predicting eccDNAs.
  -out, --output-prob   Output probability of eccDNAs present in given genomic
                        windows.
```
