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

num_elements = 5
system = Gillespie(
    num_elements,
    inits=[1, 0, 0, 0, 0],
    max_cell_num=8000
)

p0 = lambda t: (1 - 1 / (1 + np.exp(-0.6 * (t - 17))))
p1 = lambda t: (1 - 1 / (1 + np.exp(-0.6 * (t - 17))))
p2 = lambda t: (1 - 1 / (1 + np.exp(-0.6 * (t - 17))))
p3 = lambda t: (1 - 1 / (1 + np.exp(-0.6 * (t - 18))))
p4 = lambda t: (1 - 1 / (1 + np.exp(-0.6 * (t - 19))))
d0 = lambda t: 1 - p0(t)
d1 = lambda t: 1 - p1(t)
d2 = lambda t: 1 - p2(t)
d3 = lambda t: 1 - p3(t)

system.add_reaction(p0, [1, 0, 0, 0, 0], [2, 0, 0, 0, 0], index=0) # 0 self renew
system.add_reaction(p1, [0, 1, 0, 0, 0], [0, 2, 0, 0, 0], index=1) # 1 self renew
system.add_reaction(p2, [0, 0, 1, 0, 0], [0, 0, 2, 0, 0], index=2) # 2 self renew
system.add_reaction(p3, [0, 0, 0, 1, 0], [0, 0, 0, 2, 0], index=3) # 3 self renew
system.add_reaction(p4, [0, 0, 0, 0, 1], [0, 0, 0, 0, 2], index=4) # 4 self renew
system.add_reaction(d0, [1, 0, 0, 0, 0], [0, 1, 0, 0, 0], index=5) # 0 -> 1 differentiation
system.add_reaction(d1, [0, 1, 0, 0, 0], [0, 0, 1, 0, 0], index=6) # 1 -> 2 differentiation
system.add_reaction(d2, [0, 0, 1, 0, 0], [0, 0, 0, 1, 0], index=7) # 2 -> 3 differentiation
system.add_reaction(d3, [0, 0, 0, 1, 0], [0, 0, 0, 0, 1], index=8) # 3 -> 4 differentiation

system.evolute(20000000)

data_path = "../results/"
tree_file_name = f"tree_origin_linear_{filename}.csv"
cell_num_file_name = f"cell_num_linear_{filename}.csv"

curr_cells = []
t = np.array(system.generation_time)
cell_num_traj = np.array(system.n)

for i in system.curr_cells.values():
    curr_cells += i
# while tree_file_name in os.listdir(data_path):
#     tree_file_name = tree_file_name[:-1] + str(int(tree_file_name[-1]) + 1)
#     cell_num_file_name = cell_num_file_name[:-1] + str(int(cell_num_file_name[-1]) + 1)

np.savetxt(
    data_path + cell_num_file_name,
    np.hstack((t.reshape(-1, 1), cell_num_traj)),
    fmt="%.5f",
)

sim_utils.wirte_lineage_info(
    data_path + tree_file_name, system.anc_cells, curr_cells, system.t[-1]
)
try:
    reconstruct(f'../results/tree_origin_linear_{filename}.csv', output=f'../results/linear_tree_gt_{filename}.nwk', num=1000, is_balance=True)
except:
    os.system(f'rm ../results/{tree_file_name}')
    os.system(f'rm ../results/{cell_num_file_name}')

tree_file = f'../results/linear_tree_gt_{filename}.nwk'
phylo_tree, branch_colors = loadtree(tree_file)
sampled_cells = [i.name for i in phylo_tree.get_terminals()]
# cell_names, cell_states, cell_generation = get_annotation('../results/tree_origin_linear.csv')
# cell_states = pd.DataFrame(data=cell_states, index=cell_names).loc[sampled_cells]
# cell_generation = pd.DataFrame(data=cell_generation, index=cell_names).loc[sampled_cells].to_numpy()

seqs = dna_seq_simulation(phylo_tree, mut_prob=0.001)

mt_muts, mutid = mtmutation(phylo_tree, mut_rate=0.0016, init_mut_rate=5, mt_copynumber=2)

pre_existing_mut = set()
for i in mt_muts['<0_0>']:
    pre_existing_mut = pre_existing_mut.union(i)
    
mf = mut_freq(mt_muts, mutid, sel_cells=[i.name for i in phylo_tree.get_terminals()])

mut_table = mf.astype(bool).astype(int)
with open(f'../results/mt_constant_mut_{filename}.phy', 'w') as f:
    f.write('{} {}\n'.format(*mut_table.shape))
    for cell in mut_table.index:
        f.write('{} {}\n'.format(cell[1:-1], ''.join(mut_table.loc[cell].astype(str)).replace('0', 'A').replace('1', 'G')))

with open(f'../results/dna_constant_mut_{filename}.phy', 'w') as f:
    f.write('{} {}\n'.format(*seqs.shape))
    for cell in seqs.index:
        f.write('{} {}\n'.format(cell[1:-1], ''.join(seqs.loc[cell].astype(str)).replace('0', 'A').replace('1', 'G')))