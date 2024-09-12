import os
import datetime
import warnings

sl_script = [
'#!/bin/bash',
f'#SBATCH -J mt_sim30',
'#SBATCH -p all',
'#SBATCH -N 1',
'#SBATCH -n 1',
'#SBATCH --mem=10G',
'#SBATCH -t 0',
'#SBATCH -o oe/%x-%j.log ',
'#SBATCH -e oe/%x-%j.err' ]


path = '/data3/wangkun/mtsim_res/240829'

for _ in range(20):
    simid = str(datetime.datetime.now().microsecond)
    with open('./sscript', 'w') as f:
        for l in sl_script:
            f.write(f'{l}\n')
        f.write(f"python sim_30gen.py -p {path} -i {simid} ")
    os.system('sbatch sscript')
