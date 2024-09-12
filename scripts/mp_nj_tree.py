import os
import pickle
import sys
sys.path.append('..')
from mtDNAsim import *
import argparse


parser = argparse.ArgumentParser(description=' ')
parser.add_argument('-i', type=str)
parser.add_argument('-p', type=str)


path = parser.parse_args().p
simid = parser.parse_args().i
if 'const' in path:
    bn = 'const'
else:
    bn = 'mid'
mu = 0.8

for gen in [400]:
    for p in [0.1, 0.9]:
        mts = pickle.load(open(f'{path}/{simid}/mt_allmuts_{bn}_0.1_{simid}_{gen}_{mu}_{p}_rs.pkl', 'rb'))
        mts = rs_cvt(mts)
        mt_freq = sparse_freq(mts, df=True)
        gt_tree = loadtree(f'{path}/{simid}/clonal_expansion_tree_{p}_{gen}_{mu}.nwk')[0]
        for i in gt_tree.get_terminals():
            if i.name.count('_')==1:
                i.name = f'{i.name}_0'
        mt_freq = mt_freq.loc[[i.name for i in gt_tree.get_terminals()]]
        for cutoff in [0, 0.01]:
            muts = mt_freq>cutoff
            muts = muts.astype(int).astype(str)
            translation_table = str.maketrans({'1': 'A', '0': 'G'})
            seqs = f'{mt_freq.shape[0]} {mt_freq.shape[1]}\n'
            for i in range(mt_freq.shape[0]):
                seqs += f'{mt_freq.index[i]} '
                seqs += ''.join(muts.iloc[i].to_numpy()).translate(translation_table)
                seqs += '\n'
            with open(f'{path}/{simid}/mt_allmuts_{bn}_0.1_{simid}_{gen}_{mu}_{p}_rs_{cutoff}.phy', 'w') as f:
                f.write(seqs)
            fn = f'/mt_allmuts_{bn}_0.1_{simid}_{gen}_{mu}_{p}_rs_{cutoff}.phy'
            os.system(f'Rscript /home/wangkun/mtDNA_simulation/scripts/tree_reconstruct.r -p {path}/{simid} -f {fn}')
        
        