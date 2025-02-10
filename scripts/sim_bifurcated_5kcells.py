import sys
sys.path.append('../..')
from mtDNAsim import *
import pandas as pd
import numpy as np
from copy import deepcopy
from collections import Counter
import argparse
import phylovelo as pv
import pickle

parser = argparse.ArgumentParser(description='save file name')
parser.add_argument('-p', type=str, help='data path')
parser.add_argument('-f', type=str, help='filename')
parser.add_argument('-bn', type=int, help='1:with/0:without bottleneck')
parser.add_argument('-nmt', type=int, nargs='+', help='number of mtDNA')
parser.add_argument('-imr', type=int, nargs='+', default=100, help='initial mutation burden')
parser.add_argument('-mu', type=float, nargs='+', help='mtDNA mutation rate per site')

filename = parser.parse_args().f
data_path = parser.parse_args().p

bn = parser.parse_args().bn
nmt = parser.parse_args().nmt
imr = parser.parse_args().imr
mutrate = parser.parse_args().mu
if not hasattr(nmt, '__iter__'):
    nmt = [nmt]
if not hasattr(imr, '__iter__'):
    imr = [imr]    
if not hasattr(mutrate, '__iter__'):
    mutrate = [mutrate]

model = 'bif'
n_mtsite = 16569
if bn:
    bn = 'mid'
else:
    bn = 'const'

try:
    os.mkdir(f"{data_path}/{model}_{bn.replace('mid', '')}/")
except:
    pass

os.mkdir(f"{data_path}/{model}_{bn.replace('mid', '')}/{filename}/")

if data_path is None:
    data_path = "../results"
    
num_elements = 5

while True:
    num_elements = 5
    system = Gillespie(
        num_elements,
        inits=[1, 0, 0, 0, 0],
        max_cell_num=20000
    )

    p0 = lambda t: (1 - 1 / (1 + np.exp(-0.5 * (t - 19))))
    p1 = lambda t: (1 - 1 / (1 + np.exp(-0.6 * (t - 18))))
    p2 = lambda t: 1
    d00 = lambda t: 0.4* (1 / (1 + np.exp(-0.5 * (t - 16))))
    d01 = lambda t: 0.4* (1 / (1 + np.exp(-0.5 * (t - 16))))
    d10 = lambda t: 0.4* (1 / (1 + np.exp(-0.6 * (t - 16))))
    d11 = lambda t: 0.4* (1 / (1 + np.exp(-0.6 * (t - 16))))


    system.add_reaction(p0, [1, 0, 0, 0, 0], [2, 0, 0, 0, 0], index=0)
    system.add_reaction(p1, [0, 1, 0, 0, 0], [0, 2, 0, 0, 0], index=1)
    system.add_reaction(p1, [0, 0, 1, 0, 0], [0, 0, 2, 0, 0], index=2)
    system.add_reaction(p2, [0, 0, 0, 1, 0], [0, 0, 0, 2, 0], index=3)
    system.add_reaction(p2, [0, 0, 0, 0, 1], [0, 0, 0, 0, 2], index=4)
    system.add_reaction(d00, [1, 0, 0, 0, 0], [0, 1, 0, 0, 0], index=5)
    system.add_reaction(d01, [1, 0, 0, 0, 0], [0, 0, 1, 0, 0], index=6)
    system.add_reaction(d10, [0, 1, 0, 0, 0], [0, 0, 0, 1, 0], index=7)
    system.add_reaction(d11, [0, 1, 0, 0, 0], [0, 0, 0, 0, 1], index=8)

    system.evolute(20000000)
    
    tree_file_name = f"tree_origin_bif_{filename}.csv"
    cell_num_file_name = f"cell_num_bif_{filename}.csv"

    curr_cells = []
    t = np.array(system.generation_time)
    cell_num_traj = np.array(system.n)

    for i in system.curr_cells.values():
        curr_cells += i

    np.savetxt(
        f"{data_path}/{model}_{bn.replace('mid', '')}/{filename}/{cell_num_file_name}",
        np.hstack((t.reshape(-1, 1), cell_num_traj)),
        fmt="%.5f",
    )

    sim_utils.wirte_lineage_info(
        f"{data_path}/{model}_{bn.replace('mid', '')}/{filename}/{tree_file_name}", system.anc_cells, curr_cells, system.t[-1]
    )
    try:
        reconstruct(f"{data_path}/{model}_{bn.replace('mid', '')}/{filename}/tree_origin_bif_{filename}.csv", output=f"{data_path}/{model}_{bn.replace('mid', '')}/{filename}/bif_tree_gt_{filename}.nwk", num=5000, is_balance=True)
        break
    except:
        os.system(f"rm {data_path}/{model}_{bn.replace('mid', '')}/{filename}/{tree_file_name}")
        os.system(f"rm {data_path}/{model}_{bn.replace('mid', '')}/{filename}/{cell_num_file_name}")

tree_file = f"{data_path}/{model}_{bn.replace('mid', '')}/{filename}/bif_tree_gt_{filename}.nwk"
phylo_tree, branch_colors = loadtree(tree_file)
sampled_cells = [i.name for i in phylo_tree.get_terminals()]
cell_names, cell_states, cell_generation = get_annotation(f"{data_path}/{model}_{bn.replace('mid', '')}/{filename}/tree_origin_bif_{filename}.csv")
cell_states = pd.DataFrame(data=cell_states, index=cell_names).loc[sampled_cells]
cell_generation = pd.DataFrame(data=cell_generation, index=cell_names).loc[sampled_cells].to_numpy()

sd = scData(
    phylo_tree=phylo_tree,
    cell_states=cell_states.to_numpy().T[0].astype('int'),
    cell_generation=cell_generation.T[0].astype('int'),
    cell_names=sampled_cells
)



ge, base_expr = sim_base_expr(sd.phylo_tree,
                                 cell_states,
                                 Ngene=4000,
                                 r_variant_gene=0.2,
                                 diff_map={0:[0],1:[0],2:[0],3:[1],4:[1]},
                                 pseudo_state_time={0:[0,5], 1:[7,12], 2:[7,12], 3:[13,18], 4:[13,18]},
                                 forward_map={},
                                 mu0_loc=0,
                                 mu0_scale=1,
                                 drift_loc=0,
                                 drift_scale=0.05,
                                )

sd.count = get_count_from_base_expr(add_lineage_noise(sd.phylo_tree, base_expr), alpha=0.4)
sd.count.to_csv(f"{data_path}/{model}_{bn.replace('mid', '')}/{filename}/count_bif_{filename}.csv")
sd.dimensionality_reduction(method='tsne', scale=50, perplexity=100, target='count')
sd.Xdr.to_csv(f"{data_path}/{model}_{bn.replace('mid', '')}/{filename}/tsne_bif_{filename}.csv")

mt_cn = {
    'mid':lambda x: 1.52 if x <= 10 else (2.85 if x <= 20 else 2),
    'const':lambda x: 2 
}

for mu in mutrate:
    for imr1 in imr:
        for nmts in nmt:
            success = 0
            while not success:
                try:
                    mt_muts, mutid = mtmutation(phylo_tree, mut_rate=mu*n_mtsite, init_mut_rate=imr1/(nmts*2), mt_copynumber=mt_cn[bn], nmts=nmts)
                    success = 1
                except:
                    None
            pre_existing_mut = set()
            for i in mt_muts['<0_0>']:
                pre_existing_mut = pre_existing_mut.union(i)
            pickle.dump(mt_muts, open(f"{data_path}/{model}_{bn.replace('mid', '')}/{filename}/mt_allmuts_{bn}_{imr1}_{nmts}_{mu}_{filename}.pkl", 'wb'))
            curr_cells = dict()
            for i in phylo_tree.get_terminals():
                curr_cells[i.name] = mt_muts[i.name]
            mf = sparse_freq(curr_cells)
            mt_seqs = mf.astype(bool).astype(int)
            mt_pre = mt_seqs[mt_seqs.columns[np.isin(mt_seqs.columns, list(pre_existing_mut))]]
            mt_dn =  mt_seqs[mt_seqs.columns[~np.isin(mt_seqs.columns, list(pre_existing_mut))]]
            mf.to_csv(f"{data_path}/{model}_{bn.replace('mid', '')}/{filename}/mt_allmuts_{bn}_{imr1}_{nmts}_{mu}_{filename}.csv")


