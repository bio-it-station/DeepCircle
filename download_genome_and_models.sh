#!/bin/bash

hg38="https://osf.io/download/wzehu/?view_only=b72484b658874f9383cf25d209bc4000"
pretrained_DNABERT="https://osf.io/download/9sjhd/?view_only=b72484b658874f9383cf25d209bc4000"
fine_tuned_DNABERT="https://osf.io/download/xvzrn/?view_only=b72484b658874f9383cf25d209bc4000"

wget $hg38 -O hg38.zip
unzip hg38.zip
cp human/* ./DNABERT/eccdna/genome/human
rm -r human
rm hg38.zip

wget $pretrained_DNABERT -O pretrained_DNABERT.zip
unzip pretrained_DNABERT.zip
mv 6-new-12w-0 ./DNABERT/
rm pretrained_DNABERT.zip

wget $fine_tuned_DNABERT -O fine_tuned_DNABERT.zip
unzip fine_tuned_DNABERT.zip
cp -r ft/* ./DNABERT/examples/ft/ && rm -r ft
rm fine_tuned_DNABERT.zip

