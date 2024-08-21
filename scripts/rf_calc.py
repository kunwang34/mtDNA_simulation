import os
import pickle
import sys
sys.path.append('..')
from mtDNAsim import *
from ete3 import Tree
from Bio import Phylo
from io import StringIO
import re
import argparse
parser = argparse.ArgumentParser(description=' ')
parser.add_argument('-p', type=str)
parser.add_argument('-i', type=str)


path = parser.parse_args().p
simid = parser.parse_args().i

mutrate = '0.8'
def robinson_foulds(tree1:'Bio.Phylo.BaseTree', tree2:'Bio.Phylo.BaseTree'):
    f = StringIO()
    Phylo.write(tree1, f, 'newick')
    tree1 = Tree(f.getvalue(), format=1)
    # tree1.prune(tree1.get_leaf_names())
    f = StringIO()
    Phylo.write(tree2, f, 'newick')
    tree2 = Tree(f.getvalue(), format=1)
    # tree2.prune(tree2.get_leaf_names())
    try:
        return tree1.robinson_foulds(tree2)
    except:
        return tree1.robinson_foulds(tree2, unrooted_trees=True)
if 'const' in path:
    bn = 'const'
else:
    bn = 'mid'
# tree1 = loadtree(f'{path}/{t1}')[0]
# tree2 = loadtree(f'{path}/{t2}')[0]
# res = robinson_foulds(tree1, tree2)
with open(f'{path}/{simid}/rf_dist.txt', 'w') as f:
    f.write(f'dat\tmeth\tgen\ts\tcutoff\trf\tmaxrf\n')
    
with open(f'{path}/{simid}/rf_dist.txt', 'a') as f:
    for tree_meth in ['nj', 'mp', 'ml']:
        for gen in [100, 400]:
            for s in [0.1, 0.9]:
                try:
                    tree_gt = loadtree(f'{path}/{simid}/clonal_expansion_tree_{s}_{gen}.nwk')[0]
                    if tree_meth == 'ml':
                        tree_rec = loadtree(f'{path}/{simid}/dna_mut{mutrate}_{s}_{gen}_iqt.phy.treefile')[0]
                        for i in tree_rec.get_terminals():
                            ns = i.name.split('_')
                            i.name = f'<{ns[1]}_{ns[2]}>_{ns[-1]}'
                    else:
                        tree_rec = loadtree(f'{path}/{simid}/dna_mut{mutrate}_{s}_{gen}.phy_{tree_meth}.nwk')[0]
                    rf = robinson_foulds(tree_gt, tree_rec)
                    f.write(f'nDNA\t{tree_meth}\t{gen}\t{s}\t0\t{rf[0]}\t{rf[1]}\n')
                except:
                    pass
                for cf in [0, 0.01]: 
                    try:
                        if tree_meth == 'ml':
                            tree_rec = loadtree(f'{path}/{simid}/mt_allmuts_{bn}_0.1_{simid}_{gen}_0.4_{s}_rs_{cf}_iqt.phy.treefile')[0]
                            for i in tree_rec.get_terminals():
                                ns = i.name.split('_')
                                i.name = f'<{ns[1]}_{ns[2]}>_{ns[-1]}'
                        else:
                            tree_rec = loadtree(f'{path}/{simid}/mt_allmuts_{bn}_0.1_{simid}_{gen}_0.4_{s}_rs_{cf}.phy_{tree_meth}.nwk')[0]
                        rf = robinson_foulds(tree_gt, tree_rec)
                        f.write(f'mtDNA\t{tree_meth}\t{gen}\t{s}\t{cf}\t{rf[0]}\t{rf[1]}\n')
                    except:
                        pass


# path = '/data3/wangkun/mtsim_res/240705/linear_const'
# simid = os.listdir(path)
# for tree_meth in ['nj', 'mp']:
#     for gen in [100, 400]:
#         for s in [0.1, 0.9]:
#             for cf in [0, 0.01]:
#                 tree_pairs = []
#                 for i in simid:
#                     tree_gt = loadtree(f'{path}/{i}/clonal_expansion_tree_{s}_{gen}.nwk')[0]
#                     tree_rec = loadtree(f'{path}/{i}/mt_allmuts_const_0.1_{i}_{gen}_0.4_{s}_rs_{cf}.phy_{tree_meth}.nwk')[0]
#                     tree_pairs.append((tree_gt, tree_rec))
#                 with Pool(len(tree_pairs)) as p:
#                     res = list(tqdm(p.imap_unordered(robinson_foulds, tree_pairs), total=len(tree_pairs)))
#                 if gen == 100:
#                     rfs_100[f'{tree_meth}_{s}_{cf}'.replace('0.1', 'maintenance').replace('0.9', 'expansion')] = res
#                 else:
#                     rfs_400[f'{tree_meth}_{s}_{cf}'.replace('0.1', 'maintenance').replace('0.9', 'expansion')] = res