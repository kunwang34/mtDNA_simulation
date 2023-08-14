import sys
sys.path.append('..')
from mtDNAsim import *
import pandas as pd
import numpy as np
from copy import deepcopy
from collections import Counter
import argparse
parser = argparse.ArgumentParser(description='save file name')
parser.add_argument('-f', type=str)
filename = parser.parse_args().f

tree_file = f'../results/linear_tree_gt_{filename}.nwk'
phylo_tree, branch_colors = loadtree(tree_file)
sampled_cells = [i.name for i in phylo_tree.get_terminals()]
# cell_names, cell_states, cell_generation = get_annotation('../results/tree_origin_linear.csv')
# cell_states = pd.DataFrame(data=cell_states, index=cell_names).loc[sampled_cells]
# cell_generation = pd.DataFrame(data=cell_generation, index=cell_names).loc[sampled_cells].to_numpy()

seqs = DNAmutation(phylo_tree, mut_rate=0.5)
seqs = seqs.astype(int)

with open(f'../results/dna_mut0.5_{filename}.phy', 'w') as f:
    f.write('{} {}\n'.format(*seqs.shape))
    for cell in seqs.index:
        f.write('{} {}\n'.format(cell[1:-1], ''.join(seqs.loc[cell].astype(str)).replace('0', 'A').replace('1', 'G')))
        