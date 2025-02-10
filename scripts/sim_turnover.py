from copy import deepcopy
import numpy as np
import pandas as pd
import pickle
from Bio import Phylo
import argparse
import warnings
from tqdm import tqdm
import sys
sys.path.append('../..')
from mtDNAsim import *



parser = argparse.ArgumentParser(description=' ')
parser.add_argument('-p', type=str, help='data path')
parser.add_argument('-i', type=str, help='filename (simid)')
parser.add_argument('-s', type=float, help='exponsion coef')
parser.add_argument('-nmt', type=int, help='number of mtDNA')
parser.add_argument('-mu', type=float, help='mtDNA mutation rate per site')

data_path = parser.parse_args().p
simid = parser.parse_args().i
mt_mu = parser.parse_args().mu
nmt = parser.parse_args().nmt
s = parser.parse_args().s

warnings.warn(f'sim_id:{simid}', Warning)


if 'linear' in data_path:
    diff_model = 'linear'
else:
    diff_model = 'bif'

if 'const' in data_path:
    bottleneck = 'const'
else:
    bottleneck = 'mid'
    
if data_path is None:
    data_path = "./"

tree = Phylo.read(f'{data_path}/{simid}/{diff_model}_tree_gt_{simid}.nwk', format='newick')
tree_origin = pd.read_csv(f'{data_path}/{simid}/tree_origin_{diff_model}_{simid}.csv')

imb = 100
mt = pickle.load(open(f'{data_path}/{simid}/mt_allmuts_{bottleneck}_{imb}_{nmt}_{mt_mu}_{simid}.pkl', 'rb')) 
sel_cells = [i.name for i in tree.get_terminals()]
max_mut_id = max([max([max(list(i)+[0]) for i in mt[j]]+[0]) for j in sel_cells])


new_mts_1 = []
for cell in tree.get_terminals():
    new_mts_1.append({cell.name:deepcopy(mt[cell.name])})
    
gen = 0
with tqdm(total=400) as pbar:
    for i in [50, 50, 50, 50, 50, 50, 50, 50]:
        for _ in range(i):
            gen += 1
            cell_number = len(new_mts_1)
            if cell_number > 6000:
                p = 0.433
            elif cell_number < 4000:
                p = 0.6
            else:
                p = 0.5
            new_mts_1, max_mut_id = ncell_division_with_mt1(new_mts_1, max_mut_id, mt_mu, target_nmts=nmt, p=p, s=s)
            pbar.update(1)
        pickle.dump(new_mts_1, open(f'{data_path}/{simid}/mt_allmuts_{bottleneck}_{imb}_{nmt}_{mt_mu}_{simid}_{s}_{gen}.pkl', 'wb'))
        try:
            mts_new = dict()
            for i in new_mts_1:
                mts_new[list(i.keys())[0]] = list(i.values())[0]
            get_gt_tree(tree, mts_new, save_path = f'{data_path}/{simid}/mt_allmuts_{bottleneck}_{imb}_{nmt}_{mt_mu}_{simid}_{s}_{gen}.nwk', collapse=True)
        except:
            None


