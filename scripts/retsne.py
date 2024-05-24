from copy import deepcopy
import numpy as np
import pandas as pd
import pickle
from Bio import Phylo
import argparse
import warnings
from tqdm import tqdm
import phylovelo as pv

parser = argparse.ArgumentParser(description=' ')
parser.add_argument('-p', type=str)
parser.add_argument('-f', type=str)
parser.add_argument('-bn', type=str)
parser.add_argument('-m', type=str)


filename = parser.parse_args().f
data_path = parser.parse_args().p
diff_model = parser.parse_args().m
bottleneck = parser.parse_args().bn

count = pd.read_csv(f'/data3/wangkun/mtsim_res/res_0419/{data_path}{filename}/count_{diff_model}_{filename}.csv', index_col=0)
sd = pv.scData(count=count)
# sd.dimensionality_reduction(method='tsne', perplexity=30)
sd.dimensionality_reduction(method='tsne', scale=100, perplexity=150, target='count')
# sd.dimensionality_reduction(method='tsne', scale=500, perplexity=20, target='count')
sd.Xdr.index = sd.count.index
sd.Xdr.to_csv(f'/data3/wangkun/mtsim_res/res_0419/{data_path}{filename}/tsne_{diff_model}_{filename}.csv')