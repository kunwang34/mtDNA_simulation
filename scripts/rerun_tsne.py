import os

path = 'linear_const'
files = os.listdir(f"/data3/wangkun/mtsim_res/res_1113/{path.replace('mid', '')}/")
files = [i for i in files if i[0] !='b']
files = [i for i in files if i[0] !='l']
files = [i for i in files if i[0] !='r']
sl_script = [
'#!/bin/bash',
'#SBATCH -J tsne_rerun',
'#SBATCH -p all',
'#SBATCH -N 1',
'#SBATCH -n 1',
'#SBATCH --mem=20G',
'#SBATCH -t 0',
'#SBATCH -o oe/%x-%j.log ',
'#SBATCH -e oe/%x-%j.err' ]


for i in files:
    # if i in '22995,29303,50041,88021,123986,511657'.split(','):
    #     continue
    sim_res = os.listdir(f"/data3/wangkun/mtsim_res/res_1113/{path.replace('mid', '')}/{i}")
    sim_res = [i for i in sim_res if i.split('_')[-1]=='50.pkl']
    if len(sim_res) == 3:
        continue
    with open('./sscript', 'w') as f:
        for l in sl_script:
            f.write(f'{l}\n')
        f.write(f"python single_branch_simulation.py -p {path.replace('mid', '')}/ -f {i} -bn {path.split('_')[-1]} -m {path.split('_')[0].replace('bifurcated', 'bif')}")
    
    os.system('sbatch sscript')
    
