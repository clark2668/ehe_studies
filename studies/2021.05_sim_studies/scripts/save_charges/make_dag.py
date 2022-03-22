import numpy as np
import os
from pathlib import Path

output_directory_head = '/data/user/brianclark/IceCube/ehe/output/sim_cutfaraway/'

# datasets = [21218 (nue), 21220 (numu), 21221 (nutau), 21315 (muongun), 20787 (corsika)]
# datasets = [21218, 21220]
# datasets = [20787, 21218, 21220]
datasets = [20787]
for dataset in datasets:

	dag_file_name = f'dagman_{dataset}.dag'
	instructions = ""
	instructions += 'CONFIG config.dagman\n'
	instructions += "\n\n"
	with open(dag_file_name, 'w') as f:
		f.write(instructions)

	main_index = 0

	in_filelist = open(f'../list_{dataset}.txt', 'r')
	for in_file in in_filelist.readlines():

		# if main_index > 10:
		# 	continue

		in_file = in_file.rstrip('\n') # remove the newline
		in_dir, in_file = os.path.split(in_file) # split between path and filename
		in_file_noext = os.path.splitext(os.path.splitext(in_file)[0])[0] # strip twice (to get rid of .i3.zst)
		out_dir = output_directory_head + f'{dataset}/'
		out_file = f'{in_file_noext}'

		instructions = ""
		instructions += f'JOB job_{main_index} job.sub \n'
		instructions += f'VARS job_{main_index} indir="{in_dir}" infile="{in_file}" outdir="{out_dir}" outfile="{out_file}" \n\n'

		with open(dag_file_name, 'a') as f:
			f.write(instructions)

		main_index+=1

	in_filelist.close()
