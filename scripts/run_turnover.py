import os
import time
# path = 'bifurcated_mid'
# path = 'bifurcated_const'
# mode = '_random_sample'

data_dir = '/data3/wangkun/mtsim_res/241209'
# time.sleep(3600)
for folder in os.listdir(data_dir):
    path = f'{data_dir}/{folder}/'
    if 'const' in path:
        bn = 'const'
    else:
        bn = 'mid'
    sl_script = [
    '#!/bin/bash',
    '#SBATCH -J mt_400_sim',
    '#SBATCH -p all,fat',
    '#SBATCH -N 1',
    '#SBATCH -n 1',
    '#SBATCH --mem=10G',
    '#SBATCH -t 0',
    '#SBATCH -o oe/%x-%j.log ',
    '#SBATCH -e oe/%x-%j.err' ]

    for i in os.listdir(path):
        # for s in [0.1, 0.5, 0.9]:
        #     nmt = 500
        #     mut_rate = 5e-8
        #     with open('./sscript', 'w') as f:
        #         for l in sl_script:
        #             f.write(f'{l}\n')
        #         f.write(f"python sim_turnover.py -p {path} -i {i} -nmt {nmt} -mu {mut_rate} -s {s}")
        #     os.system('sbatch sscript')      
        s = 0.9
        # for s in [0.1, 0.9]:
        nmt = 500
        for mut_rate in [1e-8, 1e-7]:
            with open('./sscript', 'w') as f:
                for l in sl_script:
                    f.write(f'{l}\n')
                f.write(f"python sim_turnover.py -p {path} -i {i} -nmt {nmt} -mu {mut_rate} -s {s}")
            os.system('sbatch sscript') 
        # for s in [0.1, 0.9]:
        mut_rate = 5e-8
        for nmt in [750, 1000]:
            with open('./sscript', 'w') as f:
                for l in sl_script:
                    f.write(f'{l}\n')
                f.write(f"python sim_turnover.py -p {path} -i {i} -nmt {nmt} -mu {mut_rate} -s {s}")
            os.system('sbatch sscript')

