import os

# path = 'bifurcated_mid'
# path = 'bifurcated_const'
mode = '_random_sample'
folder = '240705'

for path in ['bifurcated_const', 'bifurcated_mid', 'linear_const', 'linear_mid']:
# for path in ['linear_mid']:
    files = os.listdir(f"/data3/wangkun/mtsim_res/{folder}/{path.replace('mid', '')}/")
    files = [i for i in files if i[0] !='b']
    files = [i for i in files if i[0] !='l']
    files = [i for i in files if i[0] !='r']
    sl_script = [
    '#!/bin/bash',
    f'#SBATCH -J mt_400_{path}_{mode}',
    '#SBATCH -p all',
    '#SBATCH -N 1',
    '#SBATCH -n 1',
    '#SBATCH --mem=10G',
    '#SBATCH -t 0',
    '#SBATCH -o oe/%x-%j.log ',
    '#SBATCH -e oe/%x-%j.err' ]


    for i in files:
        # if i in '22995,29303,50041,88021,123986,511657'.split(','):
        #     continue
        # if not i in '670220,518463,637545'.split(','):
        #     continue
        sim_res = os.listdir(f"/data3/wangkun/mtsim_res/{folder}/{path.replace('mid', '')}/{i}")
        # sim_res = [i for i in sim_res if i.split('_')[-1]=='50.pkl']
        # if len(sim_res) == 3:
        #     continue
        for s in [0.1, 0.95]:
            with open('./sscript', 'w') as f:
                for l in sl_script:
                    f.write(f'{l}\n')
                f.write(f"python single_branch_simulation{mode}.py -p {path.replace('mid', '')}/ -f {i} -bn {path.split('_')[-1]} -m {path.split('_')[0].replace('bifurcated', 'bif')} -mu 0.4 -s {s} -fo {folder}")

            os.system('sbatch sscript')

