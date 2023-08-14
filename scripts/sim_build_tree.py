import os
import datetime

fn = str(datetime.datetime.now().microsecond)
os.system(f'python sim_5k.py -f {fn}')
# os.system(f'Rscript tree_reconstruct.r -f {fn}')
# os.system(f'iqtree -s ../results/dna_mut_{fn}.phy -quiet')
# os.system(f'iqtree -s ../results/mt_mut_pre_{fn}.phy -quiet')
# os.system(f'iqtree -s ../results/mt_mut_dn_{fn}.phy -quiet')
# os.system(f'iqtree -s ../results/mt_mut_{fn}.phy -quiet')
# os.system(f'iqtree -s ../results/mt1_mut_pre_{fn}.phy -quiet')
# os.system(f'iqtree -s ../results/mt1_mut_dn_{fn}.phy -quiet')
# os.system(f'iqtree -s ../results/mt1_mut_{fn}.phy -quiet')
# os.system(f'iqtree -s ../results/mt5_mut_pre_{fn}.phy -quiet')
# os.system(f'iqtree -s ../results/mt5_mut_dn_{fn}.phy -quiet')
# os.system(f'iqtree -s ../results/mt5_mut_{fn}.phy -quiet')

# import re
# trees = []
# for i in os.listdir('../results/'):
#     if i[-6:] == 'mp.nwk':
#         trees.append(re.findall('[0-9]+', i)[-1])
# trees = set(trees)
# for fn in trees:
#     # os.system(f'python resim_dna.py -f {fn}')
#     # os.system(f'Rscript tree_reconstruct-dna.r -f {fn}')
#     os.system(f'iqtree -s ../results/dna_mut0.5_{fn}.phy -quiet')