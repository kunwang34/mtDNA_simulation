import os
import datetime
import warnings

sl_script = [
'#!/bin/bash',
f'#SBATCH -J mt_sim',
'#SBATCH -p all',
'#SBATCH -N 1',
'#SBATCH -n 1',
'#SBATCH --mem=10G',
'#SBATCH -t 0',
'#SBATCH -o oe/%x-%j.log ',
'#SBATCH -e oe/%x-%j.err' ]


path = '/data3/wangkun/mtsim_res/res_0421/'

for _ in range(20):
    # for m in ['linear', 'bifurcated']:
    #     for bn in ['const', 'mid']:
    for m in ['linear']:
        for bn in ['mid']:
            fn = str(datetime.datetime.now().microsecond)
            with open('./sscript', 'w') as f:
                for l in sl_script:
                    f.write(f'{l}\n')
                f.write(f"python sim_5k_{m.replace('bifurcated', 'bif')}.py -p {path} -f {fn} -m {m} -bn {bn}")
            os.system('sbatch sscript')
