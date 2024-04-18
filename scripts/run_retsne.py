import os

for path in ['bifurcated_const', 'bifurcated_mid']:
# for path in ['linear_const']:
    files = os.listdir(f"/data3/wangkun/mtsim_res/res_0415/{path.replace('mid', '')}/")
    files = [i for i in files if i[0] !='b']
    files = [i for i in files if i[0] !='l']
    files = [i for i in files if i[0] !='r']
    sl_script = [
    '#!/bin/bash',
    f'#SBATCH -J rerun_tsne',
    '#SBATCH -p all',
    '#SBATCH -N 1',
    '#SBATCH -n 1',
    '#SBATCH --mem=10G',
    '#SBATCH -t 0',
    '#SBATCH -o oe/%x-%j.log ',
    '#SBATCH -e oe/%x-%j.err' ]


    for i in files:
        # if not i in ['123456']:
        #     continue
        sim_res = os.listdir(f"/data3/wangkun/mtsim_res/res_0415/{path.replace('mid', '')}/{i}")

        with open('./sscript', 'w') as f:
            for l in sl_script:
                f.write(f'{l}\n')
            f.write(f"python retsne.py -p {path.replace('mid', '')}/ -f {i} -bn {path.split('_')[-1]} -m {path.split('_')[0].replace('bifurcated', 'bif')}")
        os.system('sbatch sscript')

