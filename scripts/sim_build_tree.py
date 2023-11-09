import os
import datetime

fn = str(datetime.datetime.now().microsecond)
path = '../results'
os.mkdir(f'{path}/{fn}')
path = f'{path}/{fn}'
os.system(f'python sim_5k.py -p {path} -f {fn}')
os.system(f'Rscript tree_reconstruct.r -p {path} -f {fn}')
files = os.listdir(path)
files = [i for i in files if i.endwith(f'{fn}.phy')]
for file in files:
    os.system(f'iqtree -s {path}/{file} -quiet')

