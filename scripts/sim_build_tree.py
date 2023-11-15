import os
import datetime
import warnings
fn = str(datetime.datetime.now().microsecond)
warnings.warn(f'sim_id:{fn}', Warning)
path = '/data3/wangkun/mtsim_res/res_1113/linear_const'
os.mkdir(f'{path}/{fn}')
path = f'{path}/{fn}'
os.system(f'python sim_5k.py -p {path} -f {fn}')
os.system(f'Rscript tree_reconstruct.r -p {path} -f {fn}')
files = os.listdir(path)
files = [i for i in files if i.endswith(f'{fn}.phy')]
for file in files:
    os.system(f'iqtree -s {path}/{file} -quiet')

