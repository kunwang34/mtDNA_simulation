import os

# path = 'bifurcated_mid'
# path = 'bifurcated_const'
mode = '_random_sample'

path = '/data3/wangkun/mtsim_res/240705/linear_const/'

mut_rate = 0.8
for path in ['/data3/wangkun/mtsim_res/240705/linear_const/', '/data3/wangkun/mtsim_res/240705/linear_/']:
    if 'const' in path:
        bn = 'const'
    else:
        bn = 'mid'
    sl_script = [
    '#!/bin/bash',
    '#SBATCH -J mt_400_sim',
    '#SBATCH -p all',
    '#SBATCH -N 1',
    '#SBATCH -n 1',
    '#SBATCH --mem=10G',
    '#SBATCH -t 0',
    '#SBATCH -o oe/%x-%j.log ',
    '#SBATCH -e oe/%x-%j.err' ]


    for i in os.listdir(path):
        # if i in '22995,29303,50041,88021,123986,511657'.split(','):
        #     continue
        # if not i in '670220,518463,637545'.split(','):
        #     continue
        # sim_res = os.listdir(f"{path}/{i}/mt_allmuts_0.1_{i}_{mut_rate}.pkl")
        # sim_res = os.listdir(f"/data3/wangkun/mtsim_res/{folder}/{path.replace('mid', '')}/{i}")
        for s in [0.1, 0.9]:
            with open('./sscript', 'w') as f:
                for l in sl_script:
                    f.write(f'{l}\n')
                f.write(f"python single_branch_simulation{mode}.py -p {path} -i {i} -mu {mut_rate} -s {s}")

            os.system('sbatch sscript')

