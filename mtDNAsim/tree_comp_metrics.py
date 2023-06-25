from ete3 import Tree
from Bio import Phylo
from io import StringIO
# import warnings

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
        return tree1.robinson_foulds(tree2)[0:2]
    except:
        return tree1.robinson_foulds(tree2, unrooted_trees=True)[0:2]

def to_root_length(tree, cells=None, use_branch_legnth=False):
    if cells is None:
        cells = [i.name for i in tree.get_terminals()]
    if use_branch_legnth:
        tree_depths = tree.depths()
        depths_map = dict()
        for i in tree_depths:
            depths_map[i.name] = tree_depths[i]
        return [depths_map[i[1:-1]] for i in cells]
    else:
        return [len(tree.get_path(tree.find_any(i))) for i in cells]
    
def colless_index(tree, cell_names=None):
    bals = []
    for node in tree.get_nonterminals():
        if len(node.clades) > 1:
            if len(node.clades) > 2:
                assert(f'polytomy detected: {node.name}')
            kv1, kv2 = len(node.clades[0].get_terminals()), len(node.clades[1].get_terminals())
            bals.append(np.abs(kv1-kv2))
    return np.array(bals)