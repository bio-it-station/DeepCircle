import gdown
import os

url = "https://drive.google.com/u/0/uc?id=1jiN-RlOeVvAWpEsm2PQLsSf0mRCzuRji&export=download"
output = "human.zip"
gdown.download(url, output, quiet = False)
os.system("unzip human.zip")
os.system("mv human ./DNABERT/eccdna/genome/")
os.system("rm human.zip")
