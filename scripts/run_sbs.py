import os

path = 'bif_const/'
files = os.listdir(f'/data3/wangkun/mtsim_res/res_1113/{path}')
files = [i for i in files if i[0] !='b']
files = [i for i in files if i[0] !='l']
files = [i for i in files if i[0] !='r']
sl_script = [
'#!/bin/bash',
'#SBATCH -J mt_800',
'#SBATCH -p all',
'#SBATCH -N 1',
'#SBATCH -n 1',
'#SBATCH --mem=10G',
'#SBATCH -t 0',
'#SBATCH -o oe/%x-%j.log ',
'#SBATCH -e oe/%x-%j.err' ]


for i in files:
    with open('./sscript', 'w') as f:
        for l in sl_script:
            f.write(f'{l}\n')
        f.write(f'python single_branch_simulation.py -p {path}/ -f {i} -b const -m bif')
    
    os.system('sbatch sscript')