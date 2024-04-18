from collections import defaultdict
from copy import deepcopy
from math import log
from scipy.sparse import coo_matrix
from tqdm.notebook import tqdm
from collections import Counter
import matplotlib.pyplot as plt
from tqdm.autonotebook import tqdm
import numpy as np
import pandas as pd
import scipy as sp
from scipy import stats
from scipy.special import factorial, comb
from scipy.stats import nbinom
from Bio import Phylo
from .gillespie_variable_paras import *
import warnings
from .gene_expr import *
import os
import re

def get_annotation(file):
    '''
    Get simulation data annotation 
    
    Args:
        file:
            Tree file from simulation script
    
    Return:
        list:
            Cell names
        list:
            Cell states
        list:
            Cell generations
    '''
    annotation = pd.read_csv(file)
    cell_names = []
    cell_states = []
    cell_generation = []
    for i in range(annotation.shape[0]):
        cell_generation.append(int(annotation.loc[i].generation))
        cell_names.append(
            f"<{int(annotation.loc[i].generation)}_{int(annotation.loc[i].cell_id)}>"
        )
        cell_states.append(int(annotation.loc[i].state))
    return cell_names, cell_states, cell_generation


def sim_base_expr(
    tree: "bio.phylo.tree",
    cell_states: "pd.DataFrame",
    Ngene: int,
    r_variant_gene: float,
    diff_map: dict,
    forward_map: dict = {},
    mu0_loc=20,
    mu0_scale=3,
    drift_loc=0,
    drift_scale=1,
    pseudo_state_time: dict = None,
):
    '''
    Simulation base expression
    
    Args:
        tree:
            Phylogenetic tree
        cell_states:
            DataFrame of cell types with index of cell names
        Ngene:
            Gene number
        r_variant_gene:
            Ratio of gene changes with differentiation
        diff_map:
            Differentiation relationships between different cell types
            {a:[b,c]} means 'a' is differentiated from 'b' and 'c'
        state_time:
            Pseudo time of each states
        forward_map:
            Only use in convergent model simulation
            {a:b} means 'a' will differentiated to 'b'
        mu0_loc:
            Mean of initial expression
        mu0_scale:
            Variation of initial expression
        drift_loc:
            Mean of gene drift 
        drift_scale:
            Variation of drift
    Returns:
        class:
            Gene expr program
        pd.DataFrame:
            base expression matrix
    '''
    base_expr = pd.DataFrame()
    depth = defaultdict(list)
    Nstate = len(diff_map)
    terminals_depths = tree.depths()
    for i in tree.get_terminals():
        depth[int(cell_states.loc[i.name].iloc[0])].append(int(terminals_depths[i]))

    start_time = {}
    end_time = {}
    for i in range(Nstate):
        start_time[i] = int(np.quantile(depth[i], 0.05))
        end_time[i] = int(np.quantile(depth[i], 0.95))

    t0 = start_time[0]
    for i in start_time:
        start_time[i] -= t0
        end_time[i] -= t0

    state_time = {}
    for i in range(Nstate):
        state_time[i] = [start_time[i], end_time[i]]

    pseudo_end_time = {0: end_time[0]}
    pseudo_start_time = {0: 0}
    for i in range(1, Nstate):
        pseudo_start_time[i] = (
            int(np.mean([pseudo_end_time[anc] for anc in diff_map[i]])) + 2
        )
        pseudo_end_time[i] = end_time[i] - start_time[i] + pseudo_start_time[i]

    if pseudo_state_time is None:
        pseudo_state_time = {}
        for i in range(Nstate):
            if not i in forward_map.values():
                pseudo_state_time[i] = [pseudo_start_time[i], pseudo_end_time[i]]
            else:
                pseudo_state_time[i] = [
                    pseudo_state_time[diff_map[i][0]][0] + 2,
                    pseudo_state_time[diff_map[i][0]][1] + 2,
                ]
    else:
        pseudo_start_time = {}
        pseudo_end_time = {}
        for i in pseudo_state_time:
            pseudo_start_time[i] = pseudo_state_time[i][0]
            pseudo_end_time = pseudo_state_time[i][1]

    ge = GeneExpr(
        Ngene=Ngene,
        r_variant_gene=r_variant_gene,
        diff_map=diff_map,
        forward_map=forward_map,
        state_time=pseudo_state_time,
    )

    ge.generate_genes(mu0_loc, mu0_scale, drift_loc, drift_scale)
    for cell in tree.get_terminals():
        cellstate = int(cell_states.loc[cell.name].iloc[0])
        tmp = ge.expr(
            cellstate,
            int(terminals_depths[cell])
            - t0
            - start_time[cellstate]
            + pseudo_start_time[cellstate],
        )
        base_expr = pd.concat((base_expr,pd.DataFrame(tmp,columns=[cell.name])), axis=1)
    base_expr = base_expr.T
    return ge, base_expr


def add_lineage_noise(
    tree: "bio.phylo.tree", base_expr_mat: "pd.DataFrame", scale=0.0001
):
    '''
    Simulate lineage noise
    
    Args:
        tree:
            Phylogenetic tree
        base_expr_mat:
            Base expression matrix from sim_base_expr
        scale:
            Lineage noise scale
    
    Return:
        pd.DataFrame:
            Base expression matrix with lineage noise
    
    '''
    noise = dict()
    base_expr_mat = deepcopy(base_expr_mat)
    ngene = base_expr_mat.shape[1]
    for cl in tree.get_terminals():
        path = tree.get_path(cl)
        noise[path[0].name] = np.zeros(ngene)
        for i, anc in enumerate(path):
            if not anc.name in noise:
                noise[anc.name] = np.random.normal(
                    loc=noise[path[i - 1].name], scale=scale
                )
    for i in base_expr_mat.index:
        base_expr_mat.loc[i] += noise[i]
    base_expr_mat = base_expr_mat.clip(lower=0)
    return base_expr_mat


def get_count_from_base_expr(base_expr_mat: "pd.DataFrame", alpha: int = 3):
    '''
    Draw gene expression count from base expression matrix
    
    Args:
        base_expr_mat:
            Base expression matrix
        alpha:
            Scale parameter of NB distribution
    
    Return:
        Gene count matrix
    '''
    base_expr_mat = deepcopy(base_expr_mat)
    for i in range(base_expr_mat.shape[0]):
        mu = base_expr_mat.iloc[i, :]
        mu = np.clip(mu, a_min=1e-5, a_max=None)
        sigma2 = mu + alpha * mu**2
        p, r = mu / (sigma2), mu**2 / (sigma2 - mu)
        base_expr_mat.iloc[i, :] = nbinom(r, p).rvs()
    return base_expr_mat


def get_count(paras:list):
    '''
    Draw random sample form NB distribution with paras = (r, p)
    
    Args:
        paras:
            NB parameters, [(r,p)]
    
    Return:
        int:
            Random sample
    '''
    r, p = [], []
    for para in paras:
        r.append(para[0])
        p.append(para[1] if para[1] >= 0 else 1)
    return stats.nbinom(r, p).rvs()


def reconstruct(file:str, output:str=None, seed:int=None, is_balance:bool=False, **kwargs):
    """
    Reconstruct phylogenetic tree of simulation data
    
    Args:
        file:
            Simulation file path
        output:
            Output newick file path
        seed:
            Random seed
        is_balance:
            Is all cell types' cell number equal
        ratio:
            How many cells to reconstruct
    Return:
        newick tree at output file
    """
    if seed:
        np.random.seed(seed)
    data = pd.read_csv(file)
    alives = data[data.is_alive == 1]

    if "ratio" in kwargs:
        sample_num = int(kwargs["ratio"] * alives.shape[0])
    elif "num" in kwargs:
        sample_num = kwargs["num"]
    else:
        warnings.warn("No ratio or num given, default plot all tree")
        sample_num = int(alives.shape[0])

    if is_balance:
        sample_index = np.array([])
        states = set(alives.state.to_numpy())
        for i in states:
            index = alives[alives.state == i].index.to_numpy()
            sample_index = np.concatenate(
                (
                    sample_index,
                    np.random.choice(
                        index, int(sample_num / len(states)), replace=False
                    ),
                )
            )

    else:
        index = alives.index.to_numpy()
        sample_index = np.random.choice(index, sample_num, replace=False)
    sample_index = list(sample_index)
    ##    print(data.loc[sample_index])

    new_keep = list(sample_index.copy())
    while new_keep:
        new_parents = []
        for i in new_keep:
            new_parents.append(
                data[
                    (data.generation == int(data.loc[i].generation) - 1)
                    & (data.cell_id == int(data.loc[i].parent_id))
                ].index[0]
            )
        while 0 in new_parents:
            new_parents.remove(0)
        new_keep = deepcopy(new_parents)
        sample_index.extend(new_parents)
    sample_index = np.unique(sample_index)
    sample_index = np.insert(sample_index, 0, 0)
    data = data.loc[sample_index]

    data["info"] = [
        "<{:d}_{:d}>:{}".format(
            int(data.loc[i]["generation"]),
            int(data.loc[i]["cell_id"]),
            data.loc[i]["time_death"] - data.loc[i]["time_birth"],
        )
        for i in data.index
    ]

    states = [
        "<{:d}_{:d}>:{:d}".format(
            int(data.loc[i]["generation"]),
            int(data.loc[i]["cell_id"]),
            int(data.loc[i]["state"]),
        )
        for i in data.index
    ]
    gen = data.generation.max()
    tree = []

    while gen:
        for pid in set(data[data.generation == gen].parent_id.to_numpy()):
            pair = data[
                np.all(list(zip(data.generation == gen, data.parent_id == pid)), axis=1)
            ]
            parent_index = data[
                np.all(
                    list(zip(data.generation == gen - 1, data.cell_id == pid)), axis=1
                )
            ].index[0]
            oi = data.loc[parent_index, "info"]
            if pair.shape[0] == 2:
                ni = "({}, {}){}".format(pair.iloc[0]["info"], pair.iloc[1]["info"], oi)
            else:
                ni = "({}){}".format(pair.iloc[0]["info"], oi)
            data.loc[parent_index, "info"] = ni
            data = data.drop(index=pair.index)
        gen -= 1

    with open(output, "w") as f:
        f.write(data.loc[0, "info"])
        f.write("\n")
        for i in states:
            f.write(i)
            f.write("\t")


def wirte_lineage_info(filepath, anc_cells, curr_cells, curr_time):
    '''
    Record lineage infomation in simulation
    '''
    with open(filepath, mode="w") as f:
        f.write(
            ",".join(
                [
                    "generation",
                    "cell_id",
                    "parent_id",
                    "state",
                    "is_alive",
                    "time_birth",
                    "time_death",
                ]
            )
        )
        f.write("\n")
        for c in anc_cells:
            f.write(
                ",".join(
                    [
                        str(c.gen),
                        str(c.cid),
                        str(c.parent),
                        str(c.state),
                        "0",
                        str(c.tb),
                        str(c.td),
                    ]
                )
            )
            f.write("\n")

        for c in curr_cells:
            f.write(
                ",".join(
                    [
                        str(c.gen),
                        str(c.cid),
                        str(c.parent),
                        str(c.state),
                        "1",
                        str(c.tb),
                        str(curr_time),
                    ]
                )
            )
            f.write("\n")


class Cell:
    """
    Cell class
    
    Args:
        Ngene: 
            Gene number
        state: 
            Cell type
        gen: 
            Cell generation
        cid: 
            Cell id
        parent: 
            Cell's parent
        tb: 
            Birth time
        td: 
            Death time
    """
    def __init__(
        self,
        Ngene:int = None,
        state:int = 0,
        gen:int = None,
        cid:int = None,
        parent:int = None,
        tb:float = None,
        td:float = None,
    ):
        self.state = state
        self.parent = parent
        self.cid = cid
        self.gen = gen
        self.tb = tb
        self.td = td


def dna_seq_mutate(seq, mut_prob):
    seq = np.array(seq).astype(int)
    mut_prob_arr = np.ones(len(seq))*mut_prob
    is_mut = np.random.rand(len(seq)) < mut_prob_arr
    seq[is_mut] = 1
    return seq.astype(int)
        
def dna_seq_simulation(tree, mut_prob=0.01, lseq=1e3):
    '''
    mut_prob:
        mutation probability per base per division
    lseq:
        length of DNA seq
    '''
    lseq = int(lseq)
    seqs = {tree.root.name:np.zeros(lseq)}

    for i in tqdm(Phylo.BaseTree._preorder_traverse(tree.root, lambda elem: elem.clades), total=len(tree.get_terminals())+len(tree.get_nonterminals())):
        if not i.is_terminal():
            for j in i.clades:
                seqs[j.name] = dna_seq_mutate(seqs[i.name], mut_prob)
    obs = []
    ind = []
    for i in tree.get_terminals():
        ind.append(i.name)
        obs.append(seqs[i.name])
    return pd.DataFrame(obs, index=ind)

                
def DNAmutation(tree, mut_rate=0.1):
    mutations = dict()
    global_mutid = -1
    for i in tree.get_terminals():
        mut = []
        for j in tree.get_path(i):
            if j.name in mutations:
                mut = deepcopy(mutations[j.name])
            else:
                for _ in range(np.random.poisson(mut_rate)):
                    mut.append(global_mutid+1)
                    global_mutid += 1
                mutations[j.name] = deepcopy(mut)

    mut_table = []
    cell_names = []
    for i in tree.get_terminals():
        seq = np.zeros(global_mutid+1)
        seq[mutations[i.name]]=1
        mut_table.append(seq)
        cell_names.append(i.name)
    return pd.DataFrame(np.array(mut_table), index=cell_names)

def initialize_mtmut(nmts=1000, mut_rate=0.2, len_mtdna=16569, birth_rate=1, death_rate=0.1):
    sys = Gillespie(1, [1], max_cell_num=nmts-1, mt_mut_rate=mut_rate)
    sys.add_reaction(lambda t: birth_rate, [1], [2], 0)
    sys.add_reaction(lambda t: death_rate, [1], [0], 0, rtype='death')
    sys.evolute(10000)
    return [i.mt_muts for i in sys.curr_cells[0]]

def cell_division_with_mt(mt_muts, global_mutid, mut_rate, mt_copynumber=2):
    new_mts = []
    for mt in mt_muts:
        for _ in range(int(mt_copynumber)):
            new_mts.append(deepcopy(mt))
            for new_mut in range(np.random.poisson(mut_rate)):
                new_mts[-1].add(global_mutid+1)
                global_mutid += 1
        if np.random.rand() < mt_copynumber - int(mt_copynumber):
            new_mts.append(deepcopy(mt))
            for new_mut in range(np.random.poisson(mut_rate)):
                new_mts[-1].add(global_mutid+1)
                global_mutid += 1 
    division = np.random.binomial(1, 0.5, len(new_mts)).astype(bool)
    cell1 = np.array(new_mts)[division]
    cell2 = np.array(new_mts)[~division]
    return list(cell1), list(cell2), global_mutid

def mtmutation(tree, mut_rate=1.6e-3, **params):
    '''
    Args:
        nmts (default:1000): 
            number of mts
        init_mut_rate (default:0.2): 
            initial mutation rate
        len_mtDNA (default:16569): 
            length of mtDNA seq (unable)
        birth_rate (default:1): 
            birth_rate of initial mts
        death_rate (default:0.1): 
            death_rate of initial mts
        mt_copynumber (default:2):
            copy number of mtDNA
        
    '''
    nmts = params.pop('nmts', 1000)
    init_mut_rate = params.pop('init_mut_rate', 0.2)
    len_mtDNA = params.pop('len_mtDNA', 16569)
    birth_rate = params.pop('birth_rate', 1)
    death_rate = params.pop('death_rate', 0.1)
    mt_copynumber = params.pop('mt_copynumber', 2)
    if not callable(mt_copynumber):
        mt_copynumber = lambda x: mt_copynumber
    
    # for i in range(10):
    #     try:
    init_cell_muts = initialize_mtmut(nmts, init_mut_rate, len_mtDNA, birth_rate, death_rate)
        #     break
        # except:
        #     continue
        

    global_mutid = max([max(i.union({0})) for i in init_cell_muts])
    mt_mutations = {tree.root.name:init_cell_muts}
    
    for i in tqdm(Phylo.BaseTree._preorder_traverse(tree.root, lambda elem: elem.clades), total=len(tree.get_terminals())+len(tree.get_nonterminals())):
        if not i.is_terminal():
            new_mts = cell_division_with_mt(mt_mutations[i.name], global_mutid, mut_rate, mt_copynumber(int(re.findall('(?<=<)[0-9]+', i.name)[0])))
            global_mutid = new_mts[2]
            for ind, j in enumerate(i.clades):
                mt_mutations[j.name] = new_mts[ind]
    # muts = []
    # cell_names = []
    # for i in tree.get_terminals():
    #     muts.append(mt_mutations[i.name])
    #     cell_names.append(i.name)
    return mt_mutations, global_mutid

def mut_freq(mt_muts, max_mut_id = None, sel_cells=None):
    if sel_cells is None:
        sel_cells = list(mt_muts.keys())
    if not max_mut_id:
        max_mut_id = max([max([max(list(i)+[0]) for i in mt_muts[j]]+[0]) for j in sel_cells])
    max_mut_id += 1
    mut_freqs = []
    cell_names = []
    for cell in tqdm(sel_cells):
        mut_pos = np.zeros((len(mt_muts[cell]), max_mut_id))
        for ind, mt in enumerate(mt_muts[cell]):
            mut_pos[ind][list(mt)] = 1
        mut_freqs.append(mut_pos.sum(0)/len(mt_muts[cell]))
        cell_names.append(cell)
    mf = pd.DataFrame(mut_freqs, index=cell_names)
    mf = mf[mf.columns[mf.sum()>0]]
    return mf

def rs_cvt(mts_rs):
    mts_new = dict()
    n_cells = []
    for i in mts_rs:
        n_cells.append(len(mts_rs[i]))
    living_cells = np.array(list(mts_rs.keys()))[np.array(n_cells)!=0]
    for cell in living_cells:
        for cnt, c in enumerate(mts_rs[cell]):
            mts_new[f'{cell}_{cnt}'] = c
    return mts_new

def sparse_freq(cells, df=True):
    cell_names = list(cells.keys())
    max_mut_id = max([max([max(list(i)+[0]) for i in cells[j]]+[0]) for j in cell_names])
    
    _row, _col, _data = [], [], []
    with tqdm(total=len(cell_names)) as pbar:
        for ind, cell in enumerate(cells):
            cell_muts = sum([list(i) for i in cells[cell]], [])
            nmts = len(cells[cell])
            cnt = Counter(cell_muts)
            for mut in cnt:
                _col.append(mut)
                _row.append(ind)
                _data.append(cnt[mut]/nmts)
            pbar.update(1)
    freq = coo_matrix((_data, (_row, _col))).tocsr()
    mut_id = np.arange(freq.shape[1])
    sel = np.array(freq.sum(axis=0)!=0).flatten()
    mut_id = mut_id[sel]
    freq = freq[:, sel]
    if df:
        freq = pd.DataFrame(freq.A, index=cell_names, columns=mut_id)
    return freq
