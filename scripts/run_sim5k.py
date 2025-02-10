import os
import datetime
import warnings

sl_script = [
'#!/bin/bash',
f'#SBATCH -J mt_sim',
'#SBATCH -p all,fat',
'#SBATCH -N 1',
'#SBATCH -n 1',
'#SBATCH --mem=10G',
'#SBATCH -t 0',
'#SBATCH -o oe/%x-%j.log ',
'#SBATCH -e oe/%x-%j.err' ]


path = '/data3/wangkun/mtsim_res/241209'

for _ in range(20):
    for m in ['linear', 'bifurcated']:
        for bn in [0, 1]:
            fn = str(datetime.datetime.now().microsecond)
            with open('./sscript', 'w') as f:
                for l in sl_script:
                    f.write(f'{l}\n')
                f.write(f"python sim_{m}_5kcells.py -p {path} -f {fn} -bn {bn} -nmt 500 750 1000 -mu 1e-8 5e-8 1e-7")
            os.system('sbatch sscript')
            
