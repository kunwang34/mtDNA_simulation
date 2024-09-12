import sys
sys.path.append('..')
from mtDNAsim import *
import pandas as pd
import numpy as np
from copy import deepcopy
from collections import Counter
import argparse
import pickle

parser = argparse.ArgumentParser(description=' ')
parser.add_argument('-i', type=str)
parser.add_argument('-p', type=str)
parser.add_argument('-mu', type=float)
parser.add_argument('-nmt', type=int)

data_path = parser.parse_args().p
simid = parser.parse_args().i
mut_rate = parser.parse_args().mu
n_mts = parser.parse_args().nmt

if 'const' in data_path:
    bn = 'const'
else:
    bn = 'mid'
    
tree_file = f"{data_path}/{simid}/linear_tree_gt_{simid}.nwk"
phylo_tree, branch_colors = loadtree(tree_file)

mt_cn = {
    'mid':lambda x: 1.52 if x <= 10 else (2.85 if x <= 20 else 2),
    'const':lambda x: 2 
}

imr= 0.1
success = 0
while not success:
    try:
        mt_muts, mutid = mtmutation(phylo_tree, mut_rate=mut_rate/n_mts, init_mut_rate=imr, mt_copynumber=mt_cn[bn], nmts=n_mts)
        success = 1
    except:
        pass
pickle.dump(mt_muts, open(f"{data_path}/{simid}/mt_allmuts_{bn}_{imr}_{simid}_{mut_rate}.pkl", 'wb'))

