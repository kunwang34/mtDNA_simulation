import os

path = '/data3/wangkun/mtsim_res/240705/linear_const'
# simids = os.listdir(path)
# for path in ['/data3/wangkun/mtsim_res/240705/linear_const', '/data3/wangkun/mtsim_res/240705/linear_']:
for path in ['/data3/wangkun/mtsim_res/240705/linear_const']:
    for simid in os.listdir(path):
        sl_script = [
        '#!/bin/bash',
        '#SBATCH -J tree',
        '#SBATCH -p all',
        '#SBATCH -N 1',
        '#SBATCH -n 4',
        '#SBATCH --mem=8G',
        '#SBATCH -t 0',
        '#SBATCH -o oe/%x-%j.log ',
        '#SBATCH -e oe/%x-%j.err' ]

        with open('./sscript', 'w') as f:
            for l in sl_script:
                f.write(f'{l}\n')
            f.write(f"python DNA_tree.py -p {path} -i {simid}")
            # f.write(f"python mp_nj_tree.py -p {path} -i {simid}")
            # f.write(f"python ml_tree.py -p {path} -i {simid}")

        os.system('sbatch sscript')

