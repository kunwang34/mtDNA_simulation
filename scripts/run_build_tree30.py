import os

path = '/data3/wangkun/mtsim_res/20240903'

# for simid in os.listdir(path):
sl_script = [
'#!/bin/bash',
'#SBATCH -J tree',
'#SBATCH -p all,fat',
'#SBATCH --exclude node3',
'#SBATCH -N 1',
'#SBATCH -n 1',
'#SBATCH --mem=4G',
'#SBATCH -t 0',
'#SBATCH -o oe/%x-%j.log ',
'#SBATCH -e oe/%x-%j.err' ]
for model in ['mid']:
    for nrm in [10, 100]:
        for simid in os.listdir(f'{path}/{model}_{nrm}'):
            with open('./sscript', 'w') as f:
                for l in sl_script:
                    f.write(f'{l}\n')
                f.write(f"python build_tree_30.py -p {path}/{model}_{nrm} -i {simid}")
                # f.write(f"python mp_nj_tree.py -p {path} -i {simid}")
                # f.write(f"python ml_tree.py -p {path} -i {simid}")

            os.system('sbatch sscript')

