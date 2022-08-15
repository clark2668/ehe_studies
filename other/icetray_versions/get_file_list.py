'''
Save the path to the first data file in a run to a list 
(faster to access later)
'''


from icecube.phys_services import goodrunlist
import numpy as np

good_run_lists = {
    2010: '/data/exp/IceCube/2010/filtered/level2pass2a/IC79_2010_GoodRunInfo.txt',
    2011: '/data/exp/IceCube/2011/filtered/level2pass2a/IC86_2011_GoodRunInfo.txt',
    2012: '/data/exp/IceCube/2012/filtered/level2pass2a/IC86_2012_GoodRunInfo.txt',
    2013: '/data/exp/IceCube/2013/filtered/level2pass2a/IC86_2013_GoodRunInfo.txt',
    2014: '/data/exp/IceCube/2014/filtered/level2pass2a/IC86_2014_GoodRunInfo.txt',
    2015: '/data/exp/IceCube/2015/filtered/level2pass2a/IC86_2015_GoodRunInfo.txt',
    2016: '/data/exp/IceCube/2016/filtered/level2pass2a/IC86_2016_GoodRunInfo.txt',
    2017: '/data/exp/IceCube/2017/filtered/level2/IC86_2017_GoodRunInfo.txt',
    2018: '/data/exp/IceCube/2018/filtered/level2/IC86_2018_GoodRunInfo.txt',
    2019: '/data/exp/IceCube/2019/filtered/level2/IC86_2019_GoodRunInfo.txt',
    2020: '/data/exp/IceCube/2020/filtered/level2/IC86_2020_GoodRunInfo.txt',
    2021: '/data/exp/IceCube/2021/filtered/level2/IC86_2021_GoodRunInfo.txt',
}

years = np.arange(2010, 2022)
from tqdm import tqdm

for y in years:

    file = f'file_list_{y}.txt'
    print(f"Working on {y}")

    grl = goodrunlist.GoodRunList()
    grl.load(good_run_lists[y])
    run_ids = grl.get_run_ids()
    for run_id, data in tqdm(grl.items()):
    # for run_id, data in grl.items():
        l2_files = data.get_files()
        if len(l2_files) > 0:
            first_file = l2_files[0]

            with open(file, 'a') as f:
                f.write(f'{run_id},{first_file}\n')


