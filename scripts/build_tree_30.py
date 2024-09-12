import os
import pickle
import sys
sys.path.append('..')
from mtDNAsim import *
import argparse
from Bio import Phylo

parser = argparse.ArgumentParser(description=' ')
parser.add_argument('-i', type=str)
parser.add_argument('-p', type=str)


path = parser.parse_args().p
simid = parser.parse_args().i

translation_table = str.maketrans({'A': '1', 'G': '0'})

# for simid in os.listdir(path):
print(path,simid) 
for cutoff in [0, 0.01]:
    for gen in [30, 130, 330]:
        fn = f'mt_allmuts_{gen}_{cutoff}.phy'
        os.system(f'Rscript /home/wangkun/mtDNA_simulation/scripts/tree_reconstruct.r -p {path}/{simid} -f {fn}')
        fn = f'mt_allmuts_{gen}_{cutoff}_seq.phy'
        os.system(f'Rscript /home/wangkun/mtDNA_simulation/scripts/tree_reconstruct.r -p {path}/{simid} -f {fn}')

#     with open(f'{path}/{simid}/mt_allmuts_const_0.1_{cutoff}.phy', 'r') as f2:
#         lines = f2.readlines()
#     with open(f'{path}/{simid}/mt_allmuts_const_0.1_{cutoff}_iqt.phy', 'w') as f1:
#         f1.write('\n'.join(lines).translate(translation_table))
#     os.system(f'iqtree -s {path}/{simid}/mt_allmuts_const_0.1_{cutoff}_iqt.phy -redo -mem 10G -nt 8 -st BIN -quiet')  
    

            