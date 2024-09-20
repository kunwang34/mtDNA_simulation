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
from Bio import Phylo

parser = argparse.ArgumentParser(description=' ')
parser.add_argument('-p', type=str)
# parser.add_argument('-i', type=str)



path = parser.parse_args().p
# simid = parser.parse_args().i

for gen in [30,130,330]:
    mt_muts = pickle.load(open(f'{path}/mt_allmuts_{gen}.pkl', 'rb'))
    freq = sparse_freq(mt_muts)
    if gen == 30:
        tree_gt = Phylo.read(f'{path}/gt_tree.nwk', 'newick')
        freq = freq.loc[[i.name for i in tree_gt.get_terminals()]]
    cells = list(freq.index)
    np.random.shuffle(cells)
    freq = freq.loc[cells]
    read_cnt, freq_samp = sequence_sim(freq, 50, 2.5, True, [1,2,3,4,5])

    read_cnt.to_csv(f'{path}/sequence_readcount_{gen}.csv')
    for suf, freq in zip(['_seq1', '_seq2', '_seq3', '_seq4', '_seq5'],freq_samp):
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
            with open(f'{path}/mt_allmuts_{gen}_{cutoff}{suf}.phy', 'w') as f:
                f.write(seqs)
