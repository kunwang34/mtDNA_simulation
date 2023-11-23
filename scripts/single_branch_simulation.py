from copy import deepcopy
import numpy as np
import pandas as pd
import pickle
from Bio import Phylo
import argparse

parser = argparse.ArgumentParser(description=' ')
parser.add_argument('-p', type=str)
parser.add_argument('-f', type=str)
parser.add_argument('-b', type=str)
parser.add_argument('-m', type=str)

filename = parser.parse_args().f
data_path = parser.parse_args().p
diff_model = parser.parse_args().m
bottleneck = parser.parse_args().b

if data_path is None:
    data_path = "./"

def cell_division_with_mt1(mt_muts, global_mutid, mut_rate, mt_copynumber=2):
    new_mts = []
    if mt_copynumber == 2:
        new_mts = mt_muts*2
    elif mt_copynumber > 2:
        new_mts = mt_muts*2
        n_mts = len(mt_muts)
        addi = np.random.choice(range(n_mts), int(n_mts*(mt_copynumber-2)), replace=False)
        new_mts = new_mts + list(np.array(mt_muts)[addi])
    else:
        new_mts = mt_muts
        n_mts = len(mt_muts)
        addi = np.random.choice(range(n_mts), int(n_mts*(mt_copynumber-2)), replace=False)
        new_mts = new_mts + list(np.array(mt_muts)[addi])
        
    division = np.random.binomial(1, 0.5, len(new_mts)).astype(bool)
    if np.sum(division) < 100:
        division = ~division
    cell1 = np.array(new_mts)[division]
    cell2 = np.array(new_mts)[~division]
    return list(cell1), list(cell2), global_mutid

tree = Phylo.read(f'/data3/wangkun/mtsim_res/res_1113/{data_path}{filename}/{diff_model}_tree_gt_{filename}.nwk', format='newick')
tree_origin = pd.read_csv(f'/data3/wangkun/mtsim_res/res_1113/{data_path}{filename}/tree_origin_{diff_model}_{filename}.csv')
for imb in [0.1, 1, 5]:
    mt = pickle.load(open(f'/data3/wangkun/mtsim_res/res_1113/{data_path}{filename}/mt_allmuts_{bottleneck}_{imb}_{filename}.pkl', 'rb'))

    new_mts_300 = dict()
    new_mts_600 = dict()
    new_mts_800 = dict()
    if diff_model == 'const':
        for cell in tree.get_terminals():
            new_mts_300[cell.name] = mt[cell.name]
            for _ in range(300):
                new_mts_300[cell.name] = cell_division_with_mt1(new_mts_300[cell.name],0,0)[0]
            new_mts_600[cell.name] = new_mts_300[cell.name]
            for _ in range(300):
                new_mts_600[cell.name] = cell_division_with_mt1(new_mts_600[cell.name],0,0)[0]
            new_mts_800[cell.name] = new_mts_600[cell.name]
            for _ in range(200):
                new_mts_800[cell.name] = cell_division_with_mt1(new_mts_800[cell.name],0,0)[0]
    else:
        for cell in tree.get_terminals():
            new_mts_300[cell.name] = mt[cell.name]
            for _ in range(12):
                new_mts_300[cell.name] = cell_division_with_mt1(new_mts_300[cell.name],0,0,2.2)[0]
            for _ in range(288):
                new_mts_300[cell.name] = cell_division_with_mt1(new_mts_300[cell.name],0,0)[0]
            new_mts_600[cell.name] = new_mts_300[cell.name]
            for _ in range(300):
                new_mts_600[cell.name] = cell_division_with_mt1(new_mts_600[cell.name],0,0)[0]
            new_mts_800[cell.name] = new_mts_600[cell.name]
            for _ in range(200):
                new_mts_800[cell.name] = cell_division_with_mt1(new_mts_800[cell.name],0,0)[0]

    pickle.dump(new_mts_300, open(f'/data3/wangkun/mtsim_res/res_1113/{data_path}{filename}/mt_allmuts_{bottleneck}_{imb}_{filename}_300.pkl', 'wb'))
    pickle.dump(new_mts_600, open(f'/data3/wangkun/mtsim_res/res_1113/{data_path}{filename}/mt_allmuts_{bottleneck}_{imb}_{filename}_600.pkl', 'wb'))
    pickle.dump(new_mts_800, open(f'/data3/wangkun/mtsim_res/res_1113/{data_path}{filename}/mt_allmuts_{bottleneck}_{imb}_{filename}_800.pkl', 'wb'))