from copy import deepcopy
import numpy as np
import pandas as pd
import pickle
from Bio import Phylo
import argparse
import warnings
from tqdm import tqdm

parser = argparse.ArgumentParser(description=' ')
parser.add_argument('-p', type=str)
parser.add_argument('-f', type=str)
parser.add_argument('-fo', type=str)
parser.add_argument('-bn', type=str)
parser.add_argument('-m', type=str)
parser.add_argument('-s', type=float)
parser.add_argument('-mu', type=float, default=0.4)

filename = parser.parse_args().f
folder = parser.parse_args().fo
data_path = parser.parse_args().p
diff_model = parser.parse_args().m
bottleneck = parser.parse_args().bn
mt_mu = parser.parse_args().mu
s = parser.parse_args().s
warnings.warn(f'sim_id:{filename}', Warning)

if data_path is None:
    data_path = "./"

def cell_division_with_mt1(mt_muts, global_mutid, mut_rate, mt_copynumber=2, target_nmts=500):
    new_mts = []
    nmts = len(mt_muts)
    if nmts < target_nmts*0.8:
        mt_copynumber = 2.3
    elif nmts > target_nmts*1.2:
        mt_copynumber = 1.8
    else:
        mt_copynumber = 2
    
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
        addi = np.random.choice(range(n_mts), int(n_mts*(mt_copynumber-1)), replace=False)
        new_mts = new_mts + list(np.array(mt_muts)[addi])
        
    division = np.random.binomial(1, 0.5, len(new_mts)).astype(bool)
    while (np.sum(division) < min(0.4*len(new_mts), 100)) or (np.sum(~division) < min(0.4*len(new_mts), 100)):
        division = np.random.binomial(1, 0.5, len(new_mts)).astype(bool)
    cell1 = np.array(new_mts)[division]
    cell2 = np.array(new_mts)[~division]
    for i in np.random.choice(range(len(cell1)), np.random.poisson(mut_rate)):
        global_mutid += 1
        mt_new = deepcopy(cell1[i])
        mt_new.add(global_mutid)
        cell1[i] = mt_new
        
    for i in np.random.choice(range(len(cell2)), np.random.poisson(mut_rate)):
        global_mutid += 1
        mt_new = deepcopy(cell2[i])
        mt_new.add(global_mutid)
        cell2[i] = mt_new
    return list(cell1), list(cell2), global_mutid

def ncell_division_with_mt1(mt_muts, global_mutid, mut_rate, mt_copynumber=2, target_nmts=500, p=0.5, s=1):
    res = []
    for cell in mt_muts:
        res1, res2, global_mutid = cell_division_with_mt1(cell, global_mutid, mut_rate, mt_copynumber=mt_copynumber)
        res.append(res1)
        res.append(res2)
    n_cells = len(res)
    res_new = []
    sel = [bool(int(i)) for i in ''.join(np.random.choice('10,00,11'.split(','), n_cells//2, p=[s, (1-p)*(1-s), p*(1-s)]))]
    # for i in np.where(np.random.binomial(1, p, n_cells))[0]:
    #     res_new.append(res[i])
    for ind, choice in enumerate(sel):
        if choice:
            res_new.append(res[ind])
    return res_new, global_mutid 


tree = Phylo.read(f'/data3/wangkun/mtsim_res/{folder}/{data_path}{filename}/{diff_model}_tree_gt_{filename}.nwk', format='newick')
tree_origin = pd.read_csv(f'/data3/wangkun/mtsim_res/{folder}/{data_path}{filename}/tree_origin_{diff_model}_{filename}.csv')
# for imb in [0.1, 1, 5]:
imb = 0.1

mt = pickle.load(open(f'/data3/wangkun/mtsim_res/{folder}/{data_path}{filename}/mt_allmuts_{bottleneck}_{imb}_{filename}.pkl', 'rb'))  
sel_cells = [i.name for i in tree.get_terminals()]
max_mut_id = max([max([max(list(i)+[0]) for i in mt[j]]+[0]) for j in sel_cells])

new_mts_1 = dict()
for cell in tree.get_terminals():
    new_mts_1[cell.name] = [deepcopy(mt[cell.name])]

gen = 20
with tqdm(total=380) as pbar:
    # for i in [15, 35, 50, 200, 500]:
    for i in [30, 50, 50, 50, 50, 50, 50, 50]:
        for _ in range(i):
            gen += 1
            cell_number = np.sum([len(new_mts_1[i]) for i in new_mts_1.keys()])
            if cell_number > 6000:
                p = 0.433
            elif cell_number < 4000:
                p = 0.6
            else:
                p = 0.5
            for cell in tree.get_terminals():
                tmp = ncell_division_with_mt1(new_mts_1[cell.name], max_mut_id, mt_mu, p=p, s=s)
                max_mut_id = tmp[-1]
                new_mts_1[cell.name] = tmp[0]  
            pbar.update(1)
        pickle.dump(new_mts_1, open(f'/data3/wangkun/mtsim_res/{folder}/{data_path}{filename}/mt_allmuts_{bottleneck}_{imb}_{filename}_{gen}_{mt_mu}_{s}_rs.pkl', 'wb'))
    
    
    