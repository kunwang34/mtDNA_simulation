import os
import sys
sys.path.append('../')
from mtDNAsim import *
import pandas as pd
import numpy as np
from copy import deepcopy
from collections import Counter
import pickle
import argparse
from io import StringIO

parser = argparse.ArgumentParser(description=' ')
parser.add_argument('-p', type=str)
# parser.add_argument('-i', type=str)
parser.add_argument('-nrm', type=float)
parser.add_argument('-bn', type=str)
parser.add_argument('-mtmr', type=float)
parser.add_argument('-nmts', type=int)


data_path = parser.parse_args().p
# simid = parser.parse_args().i
nrm = parser.parse_args().nrm
bn = parser.parse_args().bn
mt_mutrate = parser.parse_args().mtmr
n_mts = parser.parse_args().nmts


imr = nrm / 1000
selection = 0.6

num_elements = 1
success = 0
while not success:
    try:
        system = Gillespie(
            num_elements,
            inits=[1],
            max_cell_num=20000
        )

        p0 = lambda t: 0.8
        dr = lambda t: 0.33
        system.add_reaction(p0, [1], [2], index=0) # 0 self renew
        system.add_reaction(dr, [1], [0], index=13) # 3 -> 4 differentiation
        system.evolute(1000000)
        success = 1
    except:
        None

curr_cells = []

for i in system.curr_cells.values():
    curr_cells += i

sim_utils.wirte_lineage_info(
    f"{data_path}/lineage_info.csv", system.anc_cells, curr_cells, system.t[-1]
)

reconstruct(f"{data_path}/lineage_info.csv", output=f"{data_path}/gt_tree.nwk", num=1000, is_balance=True)

tree = loadtree(f'{data_path}/gt_tree.nwk')[0]
for i in tree.get_nonterminals():
    i.branch_length=1
for i in tree.get_terminals():
    i.branch_length=1
Phylo.write(tree, f'{data_path}/gt_tree.nwk', format='newick')


mt_cn = {
    'mid':lambda x: 1.52 if x <= 10 else (2.85 if x <= 20 else 2),
    'const':lambda x: 2 
}

success = 0
while not success:
    try:
        mt_muts, mutid = mtmutation(tree, mut_rate=mt_mutrate/n_mts, init_mut_rate=imr, mt_copynumber=mt_cn[bn], nmts=n_mts)
        n_root_muts = len(set(sum([list(i) for i in mt_muts['<0_0>']], [])))
        if nrm == 0:
            success = 1
        elif np.abs(n_root_muts-nrm)/nrm <= 0.2:
            success = 1
        else:
            pass
    except:
        pass

pickle.dump(mt_muts, open(f"{data_path}/mt_allmuts_30.pkl", 'wb'))


# mt_muts = pickle.load(open(f'{data_path}/mt_allmuts_{bn}_{imr}.pkl', 'rb'))
mt_freq = sparse_freq(mt_muts)
tree_gt = Phylo.read(f'{data_path}/gt_tree.nwk', format='newick')
mt_freq_leave = mt_freq.loc[[i.name for i in tree_gt.get_terminals()]]
mt_freq_samp = sequence_sim(mt_freq_leave, 50, n=2.5)
for cutoff in [0, 0.01]:
    for suf, freq in zip(['', '_seq'], [mt_freq_leave, mt_freq_samp]):
        muts = freq>cutoff
        muts = muts.iloc[:, np.where(muts.sum(0)>0)[0]]
        muts = muts.astype(int).astype(str)
        translation_table = str.maketrans({'1': 'A', '0': 'G'})
        seqs = f'{muts.shape[0]} {muts.shape[1]}\n'
        for i in range(muts.shape[0]):
            seqs += f'{muts.index[i]} '
            seqs += ''.join(muts.iloc[i].to_numpy()).translate(translation_table)
            seqs += '\n'
        with open(f'{data_path}/mt_allmuts_30_{cutoff}{suf}.phy', 'w') as f:
            f.write(seqs)

sel_cells = [i.name for i in tree.get_terminals()]
max_mut_id = max([max([max(list(i)+[0]) for i in mt_muts[j]]+[0]) for j in sel_cells])
new_mts_1 = dict()
for cell in tree.get_terminals():
    new_mts_1[cell.name] = [{cell.name:deepcopy(mt_muts[cell.name])}]
    
with tqdm(total=300) as pbar:
    gen = 30
    for iters in [100, 200]:
        for _ in range(iters):     
            gen += 1
            cell_number = np.sum([len(new_mts_1[i]) for i in new_mts_1.keys()])
            if cell_number > 1200:
                p = 0.433
            elif cell_number < 800:
                p = 0.6
            else:
                p = 0.5
            for cell in tree.get_terminals():
                tmp = ncell_division_with_mt1(new_mts_1[cell.name], max_mut_id, mt_mutrate, p=p, s=selection)
                max_mut_id = tmp[-1]
                new_mts_1[cell.name] = tmp[0]  
            pbar.update(1)

        new_mts_11 = rs_cvt(new_mts_1)
        get_gt_tree(tree_gt, new_mts_11, save_path = f"{data_path}/gt_tree_{gen}.nwk")
        pickle.dump(new_mts_11, open(f"{data_path}/mt_allmuts_{gen}.pkl", 'wb'))
        mt_freq = sparse_freq(new_mts_11)
        mt_freq_samp = sequence_sim(mt_freq, 50, n=2.5)
        
        for suf, freq in zip(['', '_seq'], [mt_freq, mt_freq_samp]):
            for cutoff in [0, 0.01]:
                muts = freq>cutoff
                muts = muts.iloc[:, np.where(muts.sum(0)>0)[0]]
                muts = muts.astype(int).astype(str)
                translation_table = str.maketrans({'1': 'A', '0': 'G'})
                seqs = f'{muts.shape[0]} {muts.shape[1]}\n'
                for i in range(muts.shape[0]):
                    seqs += f'{muts.index[i]} '
                    seqs += ''.join(muts.iloc[i].to_numpy()).translate(translation_table)
                    seqs += '\n'
                with open(f'{data_path}/mt_allmuts_{gen}_{cutoff}{suf}.phy', 'w') as f:
                    f.write(seqs)

    
    