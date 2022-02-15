import gdown
import os

url = "https://drive.google.com/u/0/uc?id=1jiN-RlOeVvAWpEsm2PQLsSf0mRCzuRji&export=download"
output = "genome.zip"
gdown.download(url, output, quiet = False)
os.system("unzip genome.zip")
os.system("rm genome.zip")
