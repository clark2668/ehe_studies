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
   
# now 
f = open(f"json_files/{which_set}.detector.l2files.json")
stuff = json.load(f)
l2_files = stuff['l2_files']


dag_file_name = f'dagman_{which_set}.dag'
instructions = ""
instructions += 'CONFIG config.dagman\n\n'

with open(dag_file_name, 'w') as f:
	f.write(instructions)

for i, f in enumerate(l2_files):

	instructions = ""
	instructions += f'JOB job_{i} job.sub \n'
	instructions += f'VARS job_{i} infile="{f}"\n\n'

	with open(dag_file_name, 'a') as f:
		f.write(instructions)
