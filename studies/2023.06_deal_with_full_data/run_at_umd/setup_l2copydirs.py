# python imports
import argparse
import json, os, pathlib
import numpy as np
import pathlib

parser = argparse.ArgumentParser()
parser.add_argument("-d", type=str, 
    dest="dataset",
    help="which dataset to read")

args = parser.parse_args()
which_set = args.dataset

# top_dir = '/data/i3store/users/baclark/datawarehouse_copy'
top_dir = '/data/user/brianclark/IceCube/EHE/datawarehouse_copy'

do_setup_year_folders = False
if do_setup_year_folders:
    years = range(2010,2023,1)
    for y in years:
        
        if y < 2017:
            
            # make the diretory structure to mirror madison
            pathlib.Path(
                os.path.join(top_dir,str(y),'filtered','level2pass2a')
                ).mkdir(parents=True, exist_ok=True)
        
        elif y > 2017:
            
            # make the diretory structure to mirror madison
            pathlib.Path(
                os.path.join(top_dir,str(y),'filtered','level2')
                ).mkdir(parents=True, exist_ok=True)
        
        elif y == 2017:
            
            pathlib.Path(
                os.path.join(top_dir,str(y),'filtered','level2pass2a')
                ).mkdir(parents=True, exist_ok=True)

            pathlib.Path(
                os.path.join(top_dir,str(y),'filtered','level2')
                ).mkdir(parents=True, exist_ok=True)
        

do_setup_run_folders = False
if do_setup_run_folders:
    
    # now 
    f = open(f"json_files/{which_set}.detector.folders.json")
    stuff = json.load(f)
    month_folders = stuff['month_folders']
    run_folders = stuff['run_folders']
    year_folders = stuff['year_folders']
    for m, r, y in zip(month_folders, run_folders, year_folders):

        the_run_no = int(r.split('Run00')[1])
        
        if the_run_no < 129523:
        
            print(f"Run {the_run_no}, use level2pass2a")
            pathlib.Path(
                os.path.join(top_dir,str(y),'filtered','level2pass2a',m,r)
            ).mkdir(parents=True, exist_ok=True)
        
        else:
            print(f"Run {the_run_no}, use level2")
            pathlib.Path(
                os.path.join(top_dir,str(y),'filtered','level2',m,r)
            ).mkdir(parents=True, exist_ok=True)
            