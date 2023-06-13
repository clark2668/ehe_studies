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

top_dir = '/data/user/brianclark/IceCube/EHE/data_full/level4_v2/IC86-2014-pass2'

f = open(f"json_files/{which_set}.detector.folders.json")
stuff = json.load(f)
run_numbers = stuff['run_folders']
for r in run_numbers:
    
    just_num = r[5:]
    
    pathlib.Path(
        os.path.join(top_dir, str(just_num))
    ).mkdir(parents=True, exist_ok=True)

