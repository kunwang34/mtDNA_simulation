from copy import deepcopy
import numpy as np
import pandas as pd
import pickle
from Bio import Phylo
import argparse
import warnings

parser = argparse.ArgumentParser(description=' ')
parser.add_argument('-p', type=str)
parser.add_argument('-f', type=str)
parser.add_argument('-bn', type=str)
parser.add_argument('-m', type=str)

filename = parser.parse_args().f
data_path = parser.parse_args().p
diff_model = parser.parse_args().m
bottleneck = parser.parse_args().bn
warnings.warn(f'sim_id:{filename}', Warning)

if data_path is None:
    data_path = "./"

def cell_division_with_mt1(mt_muts, mt_copynumber=2):
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
    return list(cell1), list(cell2)

def ncell_division_with_mt1(mt_muts, mt_copynumber=2):
    res = []
    for cell in mt_muts:
        res1, res2 = cell_division_with_mt1(cell, mt_copynumber=mt_copynumber)
        res.append(res1)
        res.append(res2)
    n_cells = len(res)
    res_new = []
    for i in np.where(np.random.binomial(1, 0.5, n_cells))[0]:
        res_new.append(res[i])
    return res_new

tree = Phylo.read(f'/data3/wangkun/mtsim_res/res_1113/{data_path}{filename}/{diff_model}_tree_gt_{filename}.nwk', format='newick')
tree_origin = pd.read_csv(f'/data3/wangkun/mtsim_res/res_1113/{data_path}{filename}/tree_origin_{diff_model}_{filename}.csv')
for imb in [0.1, 1, 5]:
    mt = pickle.load(open(f'/data3/wangkun/mtsim_res/res_1113/{data_path}{filename}/mt_allmuts_{bottleneck}_{imb}_{filename}.pkl', 'rb'))   
    new_mts_50 = dict()
    new_mts_100 = dict()
    new_mts_300 = dict()
    new_mts_800 = dict()
    if diff_model == 'const':
        for cell in tree.get_terminals():
            new_mts_50[cell.name] = [mt[cell.name]]
            for _ in range(50):
                new_mts_50[cell.name] = ncell_division_with_mt1(new_mts_50[cell.name])
            new_mts_100[cell.name] = new_mts_50[cell.name]

            for _ in range(50):
                new_mts_100[cell.name] = ncell_division_with_mt1(new_mts_100[cell.name])
            new_mts_300[cell.name] = new_mts_100[cell.name]

            for _ in range(200):
                new_mts_300[cell.name] = ncell_division_with_mt1(new_mts_300[cell.name])
            new_mts_800[cell.name] = new_mts_300[cell.name]

            for _ in range(500):
                new_mts_800[cell.name] = ncell_division_with_mt1(new_mts_800[cell.name])
    else:
        for cell in tree.get_terminals():
            new_mts_50[cell.name] = [mt[cell.name]]
            for _ in range(12):
                new_mts_50[cell.name] = ncell_division_with_mt1(new_mts_50[cell.name], 2.2)
            for _ in range(38):
                new_mts_50[cell.name] = ncell_division_with_mt1(new_mts_50[cell.name])
            new_mts_100[cell.name] = new_mts_50[cell.name]
            for _ in range(50):
                new_mts_100[cell.name] = ncell_division_with_mt1(new_mts_100[cell.name])
            new_mts_300[cell.name] = new_mts_100[cell.name]
            for _ in range(200):
                new_mts_300[cell.name] = ncell_division_with_mt1(new_mts_300[cell.name])
            new_mts_800[cell.name] = new_mts_300[cell.name]
            for _ in range(500):
                new_mts_800[cell.name] = ncell_division_with_mt1(new_mts_800[cell.name])

                
    pickle.dump(new_mts_50, open(f'/data3/wangkun/mtsim_res/res_1113/{data_path}{filename}/mt_allmuts_{bottleneck}_{imb}_{filename}_50_rs.pkl', 'wb'))
    pickle.dump(new_mts_100, open(f'/data3/wangkun/mtsim_res/res_1113/{data_path}{filename}/mt_allmuts_{bottleneck}_{imb}_{filename}_100_rs.pkl', 'wb'))
    pickle.dump(new_mts_300, open(f'/data3/wangkun/mtsim_res/res_1113/{data_path}{filename}/mt_allmuts_{bottleneck}_{imb}_{filename}_300_rs.pkl', 'wb'))
    pickle.dump(new_mts_800, open(f'/data3/wangkun/mtsim_res/res_1113/{data_path}{filename}/mt_allmuts_{bottleneck}_{imb}_{filename}_800_rs.pkl', 'wb'))