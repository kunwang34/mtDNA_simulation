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
parser.add_argument('-bn', type=str)
parser.add_argument('-m', type=str)
parser.add_argument('-mu', type=float, default=0.4)

filename = parser.parse_args().f
data_path = parser.parse_args().p
diff_model = parser.parse_args().m
bottleneck = parser.parse_args().bn
mt_mu = parser.parse_args().mu
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

def ncell_division_with_mt1(mt_muts, global_mutid, mut_rate, mt_copynumber=2, target_nmts=500, p=0.5):
    res = []
    for cell in mt_muts:
        res1, res2, global_mutid = cell_division_with_mt1(cell, global_mutid, mut_rate, mt_copynumber=mt_copynumber)
        res.append(res1)
        res.append(res2)
    n_cells = len(res)
    res_new = []
    for i in np.where(np.random.binomial(1, p, n_cells))[0]:
        res_new.append(res[i])
    return res_new, global_mutid 

tree = Phylo.read(f'/data3/wangkun/mtsim_res/res_1113/{data_path}{filename}/{diff_model}_tree_gt_{filename}.nwk', format='newick')
tree_origin = pd.read_csv(f'/data3/wangkun/mtsim_res/res_1113/{data_path}{filename}/tree_origin_{diff_model}_{filename}.csv')
# for imb in [0.1, 1, 5]:
imb = 0.1

mt = pickle.load(open(f'/data3/wangkun/mtsim_res/res_1113/{data_path}{filename}/mt_allmuts_{bottleneck}_{imb}_{filename}.pkl', 'rb'))  
sel_cells = [i.name for i in tree.get_terminals()]
max_mut_id = max([max([max(list(i)+[0]) for i in mt[j]]+[0]) for j in sel_cells])

new_mts_1 = dict()
for cell in tree.get_terminals():
    new_mts_1[cell.name] = [deepcopy(mt[cell.name])]

gen = 0
with tqdm(total=800) as pbar:
    for i in [15, 35, 50, 200, 500]:
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
                tmp = ncell_division_with_mt1(new_mts_1[cell.name], max_mut_id, mt_mu, p)
                max_mut_id = tmp[-1]
                new_mts_1[cell.name] = tmp[0]  
            pbar.update(1)
        pickle.dump(new_mts_1, open(f'/data3/wangkun/mtsim_res/res_1113/{data_path}{filename}/mt_allmuts_{bottleneck}_{imb}_{filename}_{gen}_{mt_mu}_rs.pkl', 'wb'))
    
    
    
# new_mts_15 = dict()
# new_mts_50 = dict()
# new_mts_100 = dict()
# new_mts_300 = dict()
# new_mts_800 = dict()
# if diff_model == 'const':
#     for cell in tree.get_terminals():
#         new_mts_15[cell.name] = [mt[cell.name]]
#         for _ in range(15):
#             tmp = ncell_division_with_mt1(new_mts_15[cell.name], max_mut_id, mt_mu)
#             new_mts_15[cell.name] = tmp[0]
#             global_mutid = tmp[-1]
            
#         new_mts_50[cell.name] = new_mts_15[cell.name]
#         for _ in range(35):
#             tmp = ncell_division_with_mt1(new_mts_50[cell.name], max_mut_id, mt_mu)
#             new_mts_50[cell.name] = tmp[0]
#             global_mutid = tmp[-1]
            
#         new_mts_100[cell.name] = new_mts_50[cell.name]
#         for _ in range(50):
#             tmp = ncell_division_with_mt1(new_mts_100[cell.name], max_mut_id, mt_mu)
#             new_mts_100[cell.name] = tmp[0]
#             global_mutid = tmp[-1]
            
#         new_mts_300[cell.name] = new_mts_100[cell.name]
#         for _ in range(200):
#             tmp = ncell_division_with_mt1(new_mts_300[cell.name], max_mut_id, mt_mu)
#             new_mts_300[cell.name] = tmp[0]
#             global_mutid = tmp[-1]
            
#         new_mts_800[cell.name] = new_mts_300[cell.name]
#         for _ in range(500):
#             tmp = ncell_division_with_mt1(new_mts_800[cell.name], max_mut_id, mt_mu)
#             new_mts_800[cell.name] = tmp[0]
#             global_mutid = tmp[-1]
# else:
#     for cell in tree.get_terminals():
#         new_mts_15[cell.name] = [mt[cell.name]]
#         for _ in range(12):
#             tmp = ncell_division_with_mt1(new_mts_15[cell.name], max_mut_id, mt_mu, 2.2)
#             new_mts_15[cell.name] = tmp[0]
#             global_mutid = tmp[-1]
#         for _ in range(3):
#             tmp = ncell_division_with_mt1(new_mts_15[cell.name], max_mut_id, mt_mu)
#             new_mts_15[cell.name] = tmp[0]
#             global_mutid = tmp[-1]
#         new_mts_50[cell.name] = deepcopy(new_mts_15[cell.name])
#         for _ in range(35):
#             tmp = ncell_division_with_mt1(new_mts_50[cell.name], max_mut_id, mt_mu)
#             new_mts_50[cell.name] = tmp[0]
#             global_mutid = tmp[-1]
            
#         new_mts_100[cell.name] = deepcopy(new_mts_50[cell.name])
#         for _ in range(50):
#             tmp = ncell_division_with_mt1(new_mts_100[cell.name], max_mut_id, mt_mu)
#             new_mts_100[cell.name] = tmp[0]
#             global_mutid = tmp[-1]
            
#         new_mts_300[cell.name] = deepcopy(new_mts_100[cell.name])
#         for _ in range(200):
#             tmp = ncell_division_with_mt1(new_mts_300[cell.name], max_mut_id, mt_mu)
#             new_mts_300[cell.name] = tmp[0]
#             global_mutid = tmp[-1]
            
#         new_mts_800[cell.name] = deepcopy(new_mts_300[cell.name])
#         for _ in range(500):
#             tmp = ncell_division_with_mt1(new_mts_800[cell.name], max_mut_id, mt_mu)
#             new_mts_800[cell.name] = tmp[0]
#             global_mutid = tmp[-1]

# pickle.dump(new_mts_15, open(f'/data3/wangkun/mtsim_res/res_1113/{data_path}{filename}/mt_allmuts_{bottleneck}_{imb}_{filename}_15_rs.pkl', 'wb'))
# pickle.dump(new_mts_50, open(f'/data3/wangkun/mtsim_res/res_1113/{data_path}{filename}/mt_allmuts_{bottleneck}_{imb}_{filename}_50_rs.pkl', 'wb'))
# pickle.dump(new_mts_100, open(f'/data3/wangkun/mtsim_res/res_1113/{data_path}{filename}/mt_allmuts_{bottleneck}_{imb}_{filename}_100_rs.pkl', 'wb'))
# pickle.dump(new_mts_300, open(f'/data3/wangkun/mtsim_res/res_1113/{data_path}{filename}/mt_allmuts_{bottleneck}_{imb}_{filename}_300_rs.pkl', 'wb'))
# pickle.dump(new_mts_800, open(f'/data3/wangkun/mtsim_res/res_1113/{data_path}{filename}/mt_allmuts_{bottleneck}_{imb}_{filename}_800_rs.pkl', 'wb'))