import os


path = '/data3/wangkun/mtsim_res/20240903/'
for model in os.listdir(path):
    if model == 'test':
        continue
    for simid in os.listdir(f'{path}/{model}'):
        sl_script = [
        '#!/bin/bash',
        '#SBATCH -J rf_calc',
        '#SBATCH -p all,fat',
        '#SBATCH --exclude node3',
        '#SBATCH -N 1',
        '#SBATCH -n 1',
        '#SBATCH --mem=1G',
        '#SBATCH -t 0',
        '#SBATCH -o oe/%x-%j.log ',
        '#SBATCH -e oe/%x-%j.err' ]

        with open('./sscript', 'w') as f:
            for l in sl_script:
                f.write(f'{l}\n')
            f.write(f"python rf_calc.py -p {path} -i {simid} -m {model}")

        os.system('sbatch sscript')

