import os
import datetime
import warnings

path = '/data3/wangkun/mtsim_res/20240903/'


sl_script = [
'#!/bin/bash',
f'#SBATCH -J mtsim',
'#SBATCH -p all',
'#SBATCH -N 1',
'#SBATCH -n 1',
'#SBATCH --mem=10G',
'#SBATCH -t 0',
'#SBATCH -o oe/%x-%j.log ',
'#SBATCH -e oe/%x-%j.err' ]


'''
parser = argparse.ArgumentParser(description=' ')
parser.add_argument('-p', type=str)
# parser.add_argument('-i', type=str)
parser.add_argument('-nrm', type=float)
parser.add_argument('-bn', type=str)
parser.add_argument('-mtmr', type=float)
parser.add_argument('-nmts', type=int)
'''

for bn in ['mid', 'const']:
    for nrm in [0, 10, 100]:
        path_sub = f'{path}/{bn}_{nrm}'
        try:
            os.mkdir(path_sub)
        except:
            pass
        
        for _ in range(20):
            simid = str(datetime.datetime.now().microsecond)
            path_sub = f'{path}/{bn}_{nrm}/{simid}'
            try:
                os.mkdir(path_sub)
            except:
                pass

            with open('./sscript', 'w') as f:
                for l in sl_script:
                    f.write(f'{l}\n')
                f.write(f"python mtsimulation_tree.py -p {path_sub} -nrm {nrm} -bn {bn} -mtmr 0.8 -nmts 500")
            os.system('sbatch sscript')
