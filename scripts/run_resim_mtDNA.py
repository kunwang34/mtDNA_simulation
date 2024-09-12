import os

for path in ['/data3/wangkun/mtsim_res/240705/linear_const', '/data3/wangkun/mtsim_res/240705/linear_']:
    for simid in os.listdir(path):
        sl_script = [
        '#!/bin/bash',
        '#SBATCH -J resim_mtDNA',
        '#SBATCH -p all',
        '#SBATCH -N 1',
        '#SBATCH -n 1',
        '#SBATCH --mem=4G',
        '#SBATCH -t 0',
        '#SBATCH -o oe/%x-%j.log ',
        '#SBATCH -e oe/%x-%j.err' ]

        with open('./sscript', 'w') as f:
            for l in sl_script:
                f.write(f'{l}\n')
            f.write(f"python resim_mtDNA.py -p {path} -i {simid} -mu 0.8 -nmt 500")

        os.system('sbatch sscript')

