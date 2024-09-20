import os


# for path in ['/data3/wangkun/mtsim_res/20240903/const_10', '/data3/wangkun/mtsim_res/20240903/const_100']:
for folder in os.listdir('/data3/wangkun/mtsim_res/20240903'):
    if folder =='test':
        continue
    path = f'/data3/wangkun/mtsim_res/20240903/{folder}'
    for simid in os.listdir(path):
        sl_script = [
        '#!/bin/bash',
        '#SBATCH -J sequence',
        '#SBATCH -p all,fat',
        '#SBATCH --exclude node3',
        '#SBATCH -N 1',
        '#SBATCH -n 1',
        '#SBATCH --mem=4G',
        '#SBATCH -t 0',
        '#SBATCH -o oe/%x-%j.log ',
        '#SBATCH -e oe/%x-%j.err' ]

        with open('./sscript', 'w') as f:
            for l in sl_script:
                f.write(f'{l}\n')
            f.write(f"python resim_seq.py -p {path}/{simid}")

        os.system('sbatch sscript')

