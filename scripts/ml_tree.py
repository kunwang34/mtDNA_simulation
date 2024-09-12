import os
import pickle
import sys
sys.path.append('..')
from mtDNAsim import *
import argparse


parser = argparse.ArgumentParser(description=' ')
parser.add_argument('-i', type=str)
parser.add_argument('-p', type=str)
mu = 0.8

path = parser.parse_args().p
simid = parser.parse_args().i
if 'const' in path:
    bn = 'const'
else:
    bn = 'mid'
# gen = 100
translation_table = str.maketrans({'A': '1', 'G': '0'})
for gen in [400]:
    for p in [0.1, 0.9]:
        for cutoff in [0, 0.01]:
            with open(f'{path}/{simid}/mt_allmuts_{bn}_0.1_{simid}_{gen}_{mu}_{p}_rs_{cutoff}.phy', 'r') as f2:
                lines = f2.readlines()
            with open(f'{path}/{simid}/mt_allmuts_{bn}_0.1_{simid}_{gen}_{mu}_{p}_rs_{cutoff}_iqt.phy', 'w') as f1:
                f1.write('\n'.join(lines).translate(translation_table))
            os.system(f'iqtree -s {path}/{simid}/mt_allmuts_{bn}_0.1_{simid}_{gen}_{mu}_{p}_rs_{cutoff}_iqt.phy -redo -mem 10G -nt 4 -st BIN -quiet')
        
        