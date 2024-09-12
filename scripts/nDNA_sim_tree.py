import os
import pickle
import sys
sys.path.append('..')
from mtDNAsim import *
import pandas as pd
import numpy as np
from Bio import Phylo
import argparse

parser = argparse.ArgumentParser(description=' ')
parser.add_argument('-i', type=str)
parser.add_argument('-p', type=str)


path = parser.parse_args().p
simid = parser.parse_args().i

def DNAmutation(tree, mut_rate=0.1):
    mutations = dict()
    global_mutid = -1
    for i in tree.get_terminals():
        mut = []
        for j in tree.get_path(i):
            if j in mutations:
                mut = deepcopy(mutations[j])
            else:
                for _ in range(np.random.poisson(mut_rate*j.branch_length)):
                    mut.append(global_mutid+1)
                    global_mutid += 1
                mutations[j] = deepcopy(mut)

    mut_table = []
    cell_names = []
    for i in tree.get_terminals():
        seq = np.zeros(global_mutid+1)
        seq[mutations[i]]=1
        mut_table.append(seq)
        cell_names.append(i.name)
    return pd.DataFrame(np.array(mut_table), index=cell_names)


for gen in ['', '_130', '_330']:
    tree = Phylo.read(f'{path}/{simid}/gt_tree{gen}.nwk', 'newick')
    ndna_mut = DNAmutation(tree, mut_rate=0.8)
    seqs = ndna_mut.astype(int)
    with open(f'{path}/{simid}/dna_mut{gen}.phy', 'w') as f:
        f.write('{} {}\n'.format(*seqs.shape))
        for cell in seqs.index:
            f.write('{} {}\n'.format(cell, ''.join(seqs.loc[cell].astype(str)).replace('1', 'A').replace('0', 'G')))
    os.system(f'Rscript /home/wangkun/mtDNA_simulation/scripts/tree_reconstruct.r -p {path}/{simid} -f dna_mut{gen}.phy')
        