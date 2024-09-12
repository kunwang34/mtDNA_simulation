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
mutrate = '0.8'
mt_mutrate = '0.8'
# for gen in [400]:
#     for p in [0.1, 0.9]:
#         fn = f'dna_mut{mutrate}_{p}_{gen}_{mt_mutrate}.phy'
#         os.system(f'Rscript /home/wangkun/mtDNA_simulation/scripts/tree_reconstruct.r -p {path}/{simid} -f {fn}')

translation_table = str.maketrans({'A': '0', 'G': '1'})
for gen in [400]:
    for p in [0.1, 0.9]:
        with open(f'{path}/{simid}/dna_mut{mutrate}_{p}_{gen}_{mt_mutrate}.phy', 'r') as f2:
            lines = f2.readlines()
        with open(f'{path}/{simid}/dna_mut{mutrate}_{p}_{gen}_{mt_mutrate}_iqt.phy', 'w') as f1:
            f1.write('\n'.join(lines).translate(translation_table))
        os.system(f'iqtree -s {path}/{simid}/dna_mut{mutrate}_{p}_{gen}_{mt_mutrate}_iqt.phy -redo -mem 5G -nt 4 -st BIN -quiet')

        