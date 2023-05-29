import os
from collections import defaultdict
from copy import deepcopy

import matplotlib.pyplot as plt
import numpy as np
from scipy.special import comb
from tqdm.autonotebook import tqdm

from .sim_utils import *


class Cell:
    """
    Cell class
    
    Args:
        Ngene: 
            Gene number
        mtdna:
            mtDNAs
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
        mt_muts:set = set(),
    ):
        self.state = state
        self.parent = parent
        self.cid = cid
        self.gen = gen
        self.tb = tb
        self.td = td
        self.mt_muts = mt_muts

# def initialize_mts(mt_mut_rate:float, nmts:int):
#     '''
#     Initialize mt DNAs
    
#     Args:
#         mt_mut_rate:
#             Mutation rate of mt DNA
#         nmts:
#             Number of mutations
#     '''
    
    

class Reaction:
    '''
    Cell division/differentiation type
    
    Args:
        rate:
            reaction rate function
        num_lefts:
            Cell numbers before reaction
        num_right:
            Cell numbers after reaction
        index:
            Reaction index
    '''
    def __init__(self, rate:callable=None, num_lefts:list=None, num_rights:list=None, index:int=None, rtype:str=None):
        self.rate = rate
        assert len(num_lefts) == len(num_rights)
        self.num_lefts = np.array(num_lefts)
        self.num_rights = np.array(num_rights)
        self.num_diff = self.num_rights - self.num_lefts
        self.index = index
        if rtype:
            self.type = rtype
        else:
            if 2 * sum(num_lefts) == sum(num_rights):
                self.type = "proliferate"
            else:
                self.type = "differentiation"

    def combine(self, n, s):
        return np.prod(comb(n, s))

    def propensity(self, n, t):
        return self.rate(t) * self.combine(n, self.num_lefts)


class Gillespie:
    '''
    Gillespie simulation
    
    Args:
        num_elements: 
            Cell type number
        inits: 
            Initial cell number
        max_cell_num: 
            Maximum cell number
    '''
    def __init__(
        self,
        num_elements:int,
        inits:list = None,
        init_cells:list = None,
        max_cell_num:int = 20000,
        mt_mut_rate:int = None
    ):

        assert num_elements > 0
        self.num_elements = num_elements
        self.reactions = []
        if inits is None:
            self.n = [np.ones(self.num_elements)]
        else:
            self.n = [np.array(inits)]
        self.anc_cells = []
        self.curr_cells = defaultdict(list)
        if init_cells:
            self.curr_cells[0] = init_cells
        else:
            self.curr_cells[0] = [
                Cell(gen=0, parent=-1, cid=i, tb=0) for i in range(int(self.n[0][0]))
            ]

        self.generation_time = [0]
        self.mt_mut_rate = mt_mut_rate
        self.max_cell_num = max_cell_num

    def add_reaction(self, rate:callable=None, num_lefts:list=None, num_rights:list=None, index:int=None, rtype:str=None):
        '''
        Add reactions to simulation
        
        Args:
            rate:
                reaction rate function
            num_lefts:
                Cell numbers before reaction
            num_right:
                Cell numbers after reaction
            index:
                Reaction index
        '''
        assert len(num_lefts) == self.num_elements
        assert len(num_rights) == self.num_elements
        self.reactions.append(Reaction(rate, num_lefts, num_rights, index, rtype))

    def evolute(self, steps:int):
        '''
        Run simulation
        
        Args:
            steps:
                How many steps to evolute before step
        '''
        self.t = [0]
        self.log = []
        cell_num_per_gen = {0: self.n[0][0]}
        cell_id_per_gen = {0: self.n[0][0]}
        self.global_mut_id = -1

        def proliferate(i):
            node = np.random.choice(self.curr_cells[i])
            l1, l2 = deepcopy(node), deepcopy(node)
            l1.parent, l2.parent = node.cid, node.cid
            l1.gen += 1
            l2.gen += 1
            l1.tb, l2.tb = self.t[-1], self.t[-1]
            try:
                index = cell_id_per_gen[l1.gen] + 1
            except:
                index = 0
                cell_id_per_gen.update({l1.gen: 0})
                cell_num_per_gen.update({l1.gen: 0})

            l1.cid, l2.cid = index, index + 1
            
            if self.mt_mut_rate:
                for new_mut in range(np.random.poisson(self.mt_mut_rate)):
                    self.global_mut_id += 1 
                    l1.mt_muts.add(self.global_mut_id)
                for new_mut in range(np.random.poisson(self.mt_mut_rate)):
                    self.global_mut_id += 1 
                    l2.mt_muts.add(self.global_mut_id)
       
            self.curr_cells[i].append(l1)
            self.curr_cells[i].append(l2)

            cell_num_per_gen[l1.gen] += 2
            cell_id_per_gen[l1.gen] += 2
            cell_num_per_gen[l1.gen - 1] -= 1

            self.curr_cells[i].remove(node)
            node.td = self.t[-1]
            self.anc_cells.append(node)
            
        def death(i=0):
            node = np.random.choice(self.curr_cells[i])
            cell_num_per_gen[node.gen] -= 1
            self.curr_cells[i].remove(node)
            self.anc_cells.append(node)
            
        def differentiate(i, j):

            node = np.random.choice(self.curr_cells[i])
            diffcell = deepcopy(node)
            self.curr_cells[i].remove(node)

            node.td = self.t[-1]
            diffcell.tb = self.t[-1]
            self.anc_cells.append(node)
            diffcell.gen += 1

            diffcell.state = j

            try:
                index = cell_id_per_gen[diffcell.gen] + 1
            except:
                index = 0
                cell_id_per_gen.update({diffcell.gen: 0})
                cell_num_per_gen.update({diffcell.gen: 0})

            diffcell.cid = index
            cell_id_per_gen[diffcell.gen] += 1
            cell_num_per_gen[diffcell.gen] += 1
            cell_num_per_gen[diffcell.gen - 1] -= 1
            diffcell.parent = node.cid
            self.curr_cells[diffcell.state].append(diffcell)

        class SwitchCase(object):
            def case_to_func(self, reaction):
                self.reaction = reaction
                method = getattr(self, reaction.type + "1")
                return method

            def proliferate1(self):
                proliferate(list(self.reaction.num_lefts).index(1))
            
            def death1(self):
                death()

            def differentiation1(self):
                differentiate(
                    list(self.reaction.num_lefts).index(1),
                    list(self.reaction.num_rights).index(1),
                )

        cls = SwitchCase()

        with tqdm(total=self.max_cell_num) as pbar:

            for _ in range(steps):
                all_cell_num = sum(cell_num_per_gen.values())
                pbar.update(all_cell_num - pbar.n)
                if all_cell_num > self.max_cell_num:
                    # print("\n maximum cell number reached")
                    break

                avg_generation = np.dot(
                    list(cell_num_per_gen.keys()), list(cell_num_per_gen.values())
                ) / sum(cell_num_per_gen.values())

                self.generation_time.append(avg_generation)
                A = np.array(
                    [
                        rec.propensity(self.n[-1], avg_generation)
                        for rec in self.reactions
                    ]
                )

                A0 = A.sum()
                A /= A0
                t0 = -np.log(np.random.random()) / A0
                self.t.append(self.t[-1] + t0)
                react = np.random.choice(self.reactions, p=A)

                self.log.append(react.index)
                self.n.append(self.n[-1] + react.num_diff)

                cls.case_to_func(react)()
