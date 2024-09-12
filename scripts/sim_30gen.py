import sys
sys.path.append('..')
from mtDNAsim import *
import pandas as pd
import numpy as np
from copy import deepcopy
from collections import Counter
import argparse
import pickle

parser = argparse.ArgumentParser(description='save file name')
parser.add_argument('-p', type=str)
parser.add_argument('-i', type=str)

data_path = parser.parse_args().p
simid = parser.parse_args().i

try:
    os.mkdir(f"{data_path}/{simid}")
except:
    pass

num_elements = 1
# while True:
success = 1
while success:
    try:
        system = Gillespie(
            num_elements,
            inits=[1],
            max_cell_num=18000
        )

        p0 = lambda t: 0.8
        dr = lambda t: 0.33
        system.add_reaction(p0, [1], [2], index=0) # 0 self renew
        system.add_reaction(dr, [1], [0], index=13) # 3 -> 4 differentiation
        system.evolute(1000000)
        success = 0
    except:
        pass


curr_cells = []

for i in system.curr_cells.values():
    curr_cells += i

sim_utils.wirte_lineage_info(
    f"{data_path}/{simid}/lineage_info.csv", system.anc_cells, curr_cells, system.t[-1]
)

reconstruct(f"{data_path}/{simid}/lineage_info.csv", output=f"{data_path}/{simid}/gt_tree.nwk", num=1000, is_balance=True)

tree = loadtree(f'{data_path}/{simid}/gt_tree.nwk')[0]

for i in tree.get_nonterminals():
    i.branch_length=1
for i in tree.get_terminals():
    i.branch_length=1
    
Phylo.write(tree, f'{data_path}/{simid}/gt_tree.nwk', format='newick')

imr = 0.1
mt_cn = {
    'mid':lambda x: 1.52 if x <= 10 else (2.85 if x <= 20 else 2),
    'const':lambda x: 2 
}
bn = 'const'
success = 0
while not success:
    try:
        mt_muts, mutid = mtmutation(tree, mut_rate=0.0016, init_mut_rate=imr, mt_copynumber=mt_cn[bn], nmts=500)
        success = 1
    except:
        None

pickle.dump(mt_muts, open(f"{data_path}/{simid}/mt_allmuts_{bn}.pkl", 'wb'))

mt_freq = sparse_freq(mt_muts)
for cutoff in [0, 0.01]:
    muts = mt_freq>cutoff
    muts = muts.iloc[:, np.where(muts.sum(0)>0)[0]]
    muts = muts.astype(int).astype(str)
    translation_table = str.maketrans({'1': 'A', '0': 'G'})
    seqs = f'{mt_freq.shape[0]} {mt_freq.shape[1]}\n'
    for i in range(mt_freq.shape[0]):
        seqs += f'{mt_freq.index[i]} '
        seqs += ''.join(muts.iloc[i].to_numpy()).translate(translation_table)
        seqs += '\n'
    with open(f'{data_path}/{simid}/mt_allmuts_{bn}_0.1_{cutoff}.phy', 'w') as f:
        f.write(seqs)

        
seqs = DNAmutation(tree, mut_rate=0.8)
seqs = seqs.astype(int)

with open(f'{data_path}/{simid}/ndna_0.8.phy', 'w') as f:
    f.write('{} {}\n'.format(*seqs.shape))
    for cell in seqs.index:
        f.write('{} {}\n'.format(cell, ''.join(seqs.loc[cell].astype(str)).replace('0', 'A').replace('1', 'G')))
        