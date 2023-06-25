import os
import datetime

fn = str(datetime.datetime.now().microsecond)
os.system(f'python sim1.py -f {fn}')
os.system(f'Rscript tree_reconstruct.r -f {fn}')