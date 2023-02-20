
 <img src="DeepCircle_logo.gif" height="218" width="579" style="align:center" />


---------------------------------------
Predict eccDNAs based on CNN and DNABERT
---------------------------------------
## Install

```sh
git clone --recursive https://github.com/KaiLiChang/DeepCircle.git
cd DeepCircle
conda env create -f environment.yml
conda activate DeepCircle
```

## Download and unzip reference genome, pre-trained and fine-tuned DNABERT

```
bash download_genome_and_models.sh
```

### If you want to infer eccDNA-related motifs, please download and install the MEME suite (see https://meme-suite.org/meme/doc/install.html?man_type=web)


## Quick start

## CNN
---
Train a CNN based on provided eccDNAs
```Python
python3 ./bin/DeepCircle.py CNN-train -bed ./DNABERT/eccdna/db/human/PC-3_circleseq_eccdna_filt_uniq.bed
```

Predict eccDNAs based on given bed files using pre-trained CNN models
```Python
python3 ./bin/DeepCircle.py CNN-predict -bed ./DNABERT/eccdna/db/human/PC-3_circleseq_eccdna_filt_uniq.bed -model PC-3
```

Infer eccDNA-related motifs based on prediction results of CNN
```Python
python3 ./bin/DeepCircle.py infer-motif -method CNN -data PC-3 -model PC-3
```

## DNABERT
---
Fine-tune a DNABERT model based on provided eccDNAs
```Python
python3 ./bin/DeepCircle.py DNABERT-train -bed PC-3_circleseq_eccdna_filt_uniq.bed
```

Predict eccDNAs based on given bed files using fine-tuned DNABERT models
```Python
python3 ./bin/DeepCircle.py DNABERT-predict -bed PC-3_circleseq_eccdna_filt_uniq.bed -model PC-3
```

Infer eccDNA-related motifs based on prediction results of DNABERT
```Python
python3 ./bin/DeepCircle.py infer-motif -method DNABERT -data PC-3 -model PC-3
```

## Advanced usage
---
```
usage: DeepCircle [-h]
                  {CNN-train,CNN-predict,DNABERT-train,DNABERT-predict,infer-motif}
                  ...

Wrapper function for training and predicting using CNN and DNABERT within DeepCircle

positional arguments:
  {CNN-train,CNN-predict,DNABERT-train,DNABERT-predict,infer-motif}

optional arguments:
  -h, --help            show this help message and exit
```

### Train a CNN model
```
usage: DeepCircle CNN-train [-h] [-c NUMBER_OF_CPUS] [-fa GENOME_FASTA]
                            [-gap GENOME_GAP] [-g GENOME] -bed ECCDNA_BED
                            [-r RATIO] [-epo EPOCHS] [-batch BATCH_SIZE]
                            [-dir OUTPUT_PREDICTION]

optional arguments:
  -h, --help            show this help message and exit
  -c NUMBER_OF_CPUS, --number-of-cpus NUMBER_OF_CPUS
                        Number of cpus [default: 8]
  -fa GENOME_FASTA, --genome-fasta GENOME_FASTA
                        Directory to the fasta file of hg38 genome [default:
                        ./DNABERT/eccdna/genome/human/hg38_noalt.fa]
  -gap GENOME_GAP, --genome-gap GENOME_GAP
                        Directory to the bed file of gaps of hg38 genome
                        [default:
                        ./DNABERT/eccdna/genome/human/hg38_noalt_gap.bed]
  -g GENOME, --genome GENOME
                        Directory to the file containing lengths of each
                        chromosome in hg38 genome [default:
                        ./DNABERT/eccdna/genome/human/hg38_noalt.fa.genome]
  -bed ECCDNA_BED, --eccDNA-bed ECCDNA_BED
                        Directory to the bed file containing eccDNAs [default:
                        ./DNABERT/eccdna/db/human/PC-3_circleseq_eccdna_filt_u
                        niq.bed]
  -r RATIO, --ratio RATIO
                        Ratio of data used for testing [default: 0.2]
  -epo EPOCHS, --epochs EPOCHS
                        Number of epochs for training [default: 5]
  -batch BATCH_SIZE, --batch-size BATCH_SIZE
                        Number of batch size for training [default: 32]
  -dir OUTPUT_PREDICTION, --output-prediction OUTPUT_PREDICTION
                        Directory to output prediction results of the model
                        [default: ./output/CNN]
```

### Predict eccDNAs using a pre-trained CNN model
```
usage: DeepCircle CNN-predict [-h] [-c NUMBER_OF_CPUS] [-fa GENOME_FASTA]
                              [-gap GENOME_GAP] [-g GENOME] -bed ECCDNA_BED
                              [-dir OUTPUT_PREDICTION] -model
                              {C4-2,ES2,HeLaS3,LnCap,OVCAR8,PC-3,U937,leukocytes,muscle,pool_LCN,pool_LCT}

optional arguments:
  -h, --help            show this help message and exit
  -c NUMBER_OF_CPUS, --number-of-cpus NUMBER_OF_CPUS
                        Number of cpus [default: 8]
  -fa GENOME_FASTA, --genome-fasta GENOME_FASTA
                        Directory to the fasta file of hg38 genome [default:
                        ./DNABERT/eccdna/genome/human/hg38_noalt.fa]
  -gap GENOME_GAP, --genome-gap GENOME_GAP
                        Directory to the bed file of gaps of hg38 genome
                        [default:
                        ./DNABERT/eccdna/genome/human/hg38_noalt_gap.bed]
  -g GENOME, --genome GENOME
                        Directory to the file containing lengths of each
                        chromosome in hg38 genome [default:
                        ./DNABERT/eccdna/genome/human/hg38_noalt.fa.genome]
  -bed ECCDNA_BED, --eccDNA-bed ECCDNA_BED
                        Directory to the bed file containing eccDNAs [default:
                        ./DNABERT/eccdna/db/human/PC-3_circleseq_eccdna_filt_u
                        niq.bed]
  -dir OUTPUT_PREDICTION, --output-prediction OUTPUT_PREDICTION
                        Directory to output prediction results of the model
                        [default: ./output/CNN]
  -model {C4-2,ES2,HeLaS3,LnCap,OVCAR8,PC-3,U937,leukocytes,muscle,pool_LCN,pool_LCT}, --chosen-model {C4-2,ES2,HeLaS3,LnCap,OVCAR8,PC-3,U937,leukocytes,muscle,pool_LCN,pool_LCT}
                        Choose pre-trained models for predicting eccDNAs
                        [default: ./pre_trained_models/]
```
### Fine-tune a DNABERT model
```
usage: DeepCircle DNABERT-train [-h] -bed ECCDNA_BED [-sp SPECIES] [-g GENOME]
                                [-fa GENOME_FASTA] [-gap GENOME_GAP] [-k KMER]
                                [-model_path MODEL_PATH] [-batch BATCH_SIZE]
                                [-epo EPOCHS] [-c NUMBER_OF_CPUS]

optional arguments:
  -h, --help            show this help message and exit
  -bed ECCDNA_BED, --eccDNA-bed ECCDNA_BED
                        File name of the bed file containing eccDNAs, put this
                        file in DeepCircle/DNABERT/eccdna/db/{species}/
  -sp SPECIES, --species SPECIES
                        Species from which input data derived, [default:
                        human]
  -g GENOME, --genome GENOME
                        Name of the file containing lengths of each chromosome
                        in hg38 genome, put this file in
                        DeepCircle/DNABERT/eccdna/genome/{species}/ [default:
                        hg38_noalt.fa.genome]
  -fa GENOME_FASTA, --genome-fasta GENOME_FASTA
                        Name of the fasta file of the reference genome, put
                        this file in
                        DeepCircle/DNABERT/eccdna/genome/{species}/ [default:
                        hg38_noalt_gap.fa]
  -gap GENOME_GAP, --genome-gap GENOME_GAP
                        Name of the bed file containing gaps of the reference
                        genome, put this file in
                        DeepCircle/DNABERT/eccdna/genome/{species}/ [default:
                        hg38_noalt_gap.bed]
  -k KMER, --kmer KMER  Specify k-mer used to train DNABERT [default: 6]
  -model_path MODEL_PATH, --model_path MODEL_PATH
                        Directory to the pre-trained DNABERT, [default:
                        ./DNABERT/6-new-12w-0]
  -batch BATCH_SIZE, --batch-size BATCH_SIZE
                        Number of batch size for training
  -epo EPOCHS, --epochs EPOCHS
                        Number of epochs for training
  -c NUMBER_OF_CPUS, --number-of-cpus NUMBER_OF_CPUS
                        Number of cpus
```
### Predict eccDNAs using a fine-tuned DNABERT model
```
usage: DeepCircle DNABERT-predict [-h] -bed ECCDNA_BED [-sp SPECIES]
                                  [-g GENOME] [-fa GENOME_FASTA]
                                  [-gap GENOME_GAP] [-k KMER]
                                  [-batch BATCH_SIZE] [-epo EPOCHS]
                                  [-c NUMBER_OF_CPUS]
                                  [-model {C4-2,ES2,HeLaS3,LnCap,OVCAR8,PC-3,U937,leukocytes,muscle,pool_LCN,pool_LCT}]

optional arguments:
  -h, --help            show this help message and exit
  -bed ECCDNA_BED, --eccDNA-bed ECCDNA_BED
                        File name of the bed file containing eccDNAs, put this
                        file in DeepCircle/DNABERT/eccdna/db/{species}/
  -sp SPECIES, --species SPECIES
                        Species from which input data derived, [default:
                        human]
  -g GENOME, --genome GENOME
                        Name of the file containing lengths of each chromosome
                        in hg38 genome, put this file in
                        DeepCircle/DNABERT/eccdna/genome/{species}/ [default:
                        hg38_noalt.fa.genome]
  -fa GENOME_FASTA, --genome-fasta GENOME_FASTA
                        Name of the fasta file of the reference genome, put
                        this file in
                        DeepCircle/DNABERT/eccdna/genome/{species}/ [default:
                        hg38_noalt_gap.fa]
  -gap GENOME_GAP, --genome-gap GENOME_GAP
                        Name of the bed file containing gaps of the reference
                        genome, put this file in
                        DeepCircle/DNABERT/eccdna/genome/{species}/ [default:
                        hg38_noalt_gap.bed]
  -k KMER, --kmer KMER  Specify k-mer used to train DNABERT [default: 6]
  -batch BATCH_SIZE, --batch-size BATCH_SIZE
                        Number of batch size for training
  -epo EPOCHS, --epochs EPOCHS
                        Number of epochs for training
  -c NUMBER_OF_CPUS, --number-of-cpus NUMBER_OF_CPUS
                        Number of cpus
  -model {C4-2,ES2,HeLaS3,LnCap,OVCAR8,PC-3,U937,leukocytes,muscle,pool_LCN,pool_LCT}, --chosen-model {C4-2,ES2,HeLaS3,LnCap,OVCAR8,PC-3,U937,leukocytes,muscle,pool_LCN,pool_LCT}
                        Choose pre-trained models for predicting eccDNAs
                        [default: DeepCircle/DNABERT/pre_trained_models]
```
### Infer eccDNA-related motifs using prediction results
```
usage: DeepCircle infer-motif [-h] [-method {CNN,DNABERT}] -data
                              {C4-2,ES2,HeLaS3,LnCap,OVCAR8,PC-3,U937,leukocytes,muscle,pool_LCN,pool_LCT}
                              [-n NUM_MOTIF]
                              [-model {C4-2,ES2,HeLaS3,LnCap,OVCAR8,PC-3,U937,leukocytes,muscle,pool_LCN,pool_LCT}]

optional arguments:
  -h, --help            show this help message and exit
  -method {CNN,DNABERT}
                        Choose type of model to that was used to predict
                        eccDNA [default: CNN]
  -data {C4-2,ES2,HeLaS3,LnCap,OVCAR8,PC-3,U937,leukocytes,muscle,pool_LCN,pool_LCT}, --data-type {C4-2,ES2,HeLaS3,LnCap,OVCAR8,PC-3,U937,leukocytes,muscle,pool_LCN,pool_LCT}
                        Choose type of data that was used to predict eccDNA
  -n NUM_MOTIF, --num-motif NUM_MOTIF
                        Type number of motifs to infer [default: 10]
  -model {C4-2,ES2,HeLaS3,LnCap,OVCAR8,PC-3,U937,leukocytes,muscle,pool_LCN,pool_LCT}, --chosen-model {C4-2,ES2,HeLaS3,LnCap,OVCAR8,PC-3,U937,leukocytes,muscle,pool_LCN,pool_LCT}
                        Choose pre-trained models for predicting eccDNAs
```