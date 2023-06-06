# python imports
import argparse
import json, os
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("-d", type=str, 
    dest="dataset",
    help="which dataset to read")

args = parser.parse_args()
which_set = args.dataset

top_path = '/home/brianclark/IceCube/ehe/ehe_software/ehe_code/datasets/'

f = open(os.path.join(top_path, which_set))
stuff = json.load(f)
grl_path = stuff['grl_path']

from icecube.phys_services import goodrunlist

grl = goodrunlist.GoodRunList()
grl.load(grl_path)
run_ids = grl.get_run_ids()

def parse_get_folders(gcd_file):
    
    split = os.path.split(gcd_file)
    # path is in the [0] entry, the file is in the second [1] entry
    # need to dissect the former more
    # it should look like this: /data/exp/IceCube/2010/filtered/level2pass2a/0614/Run00116048
    split = split[0].split('/')
    run_folder = split[-1]
    month_folder = split[-2]
    year_folder = split[-5]
    return month_folder, run_folder, year_folder

month_folders = []
run_folders = []
year_folders = []
the_l2_files = []

for run_id in run_ids:
    gcd_file = grl[run_id].get_gcd_file()
    print(gcd_file)

    # work out the directory structure that we need
    month_folder, run_folder, year_folder = parse_get_folders(gcd_file)
    month_folders.append(month_folder)
    run_folders.append(run_folder)
    year_folders.append(year_folder)
    
    l2_files = grl[run_id].get_files()
    for f in l2_files:
        the_l2_files.append(f)

dictionary = {
    'l2_files': the_l2_files
}
json_object = json.dumps(dictionary, indent=4)
with open(f"{which_set}.l2files.json", "w") as outfile:
    outfile.write(json_object)
del json_object, dictionary

dictionary = {
    'month_folders': month_folders,
    'run_folders': run_folders,
    'year_folders': year_folders
}
json_object = json.dumps(dictionary, indent=4)

with open(f"{which_set}.folders.json", "w") as outfile:
    outfile.write(json_object)