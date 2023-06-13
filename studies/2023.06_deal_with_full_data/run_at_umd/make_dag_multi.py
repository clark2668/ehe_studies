# python imports
import argparse
import json, os, pathlib
import numpy as np
import pathlib

parser = argparse.ArgumentParser()
parser.add_argument("-d", type=str, 
    dest="dataset", required=True,
    help="which dataset to read")

args = parser.parse_args()
which_set = args.dataset

which_to_move = 'whole_runs'
if which_to_move == 'l2files':
    
    f = open(f"json_files/{which_set}.detector.l2files.json")
    stuff = json.load(f)
    files_to_move = stuff['l2_files']
elif which_to_move == 'gcdfiles':
    
    f = open(f"json_files/{which_set}.detector.folders.json")
    stuff = json.load(f)
    files_to_move = stuff['gcd_files']

elif which_to_move == 'whole_runs':
    
    f = open(f"json_files/{which_set}.detector.folders.json")
    stuff = json.load(f)
    files_to_move = stuff['run_folder_paths']


dag_file_name = f'dagman_{which_set}.dag'
instructions = ""
instructions += 'CONFIG config.dagman\n\n'

with open(dag_file_name, 'w') as f:
    f.write(instructions)

for i, f in enumerate(files_to_move):

    instructions = ""
    instructions += f'JOB job_{i} job.sub \n'
    instructions += f'VARS job_{i} infile="{f}"\n\n'

    with open(dag_file_name, 'a') as f:
        f.write(instructions)
