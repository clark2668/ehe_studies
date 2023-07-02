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
    
f = open(f"../json_files/{which_set}.detector.folders.json")
stuff = json.load(f)
files_to_move = stuff['gcd_files']
run_numbers = stuff['run_folders']

dag_file_name = f'dagman_{which_set}.dag'
instructions = ""
instructions += 'CONFIG config.dagman\n\n'

with open(dag_file_name, 'w') as f:
    f.write(instructions)

for i, f in enumerate(files_to_move):

    run_no = run_numbers[i].split('Run00')[1]    
    instructions = ""
    instructions += f'JOB job_{i} job.sub \n'
    instructions += f'VARS job_{i} set="{which_set}" runno="{run_no}" infile="{f}"\n\n'

    with open(dag_file_name, 'a') as f:
        f.write(instructions)
