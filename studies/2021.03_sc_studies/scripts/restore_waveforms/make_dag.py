import numpy as np
import os
from pathlib import Path

output_directory_head = '/home/brian/IceCube/sc_reprocess/v2/restore_waveforms/'
exclude_atwd = True


dag_file_name = 'dagman_{}.dag'.format(exclude_atwd)
instructions = ""
instructions += 'CONFIG config.dagman\n'
instructions += "\n\n"
with open(dag_file_name, 'w') as f:
	f.write(instructions)

main_index = 0

in_filelist = open(f'../list_2015_sc2_full.txt', 'r')
for in_file in in_filelist.readlines():

	# if main_index > 3:
	# 	continue

	in_file = in_file.rstrip('\n') # remove the newline
	in_dir, in_file = os.path.split(in_file) # split between path and filename
	in_file_noext = os.path.splitext(os.path.splitext(in_file)[0])[0] # strip twice (to get rid of .i3.zst/i3.bz2)
	out_dir = output_directory_head
	out_file = f'{in_file_noext}'

	instructions = ""
	instructions += f'JOB job_{main_index} job.sub \n'
	instructions += f'VARS job_{main_index} indir="{in_dir}" infile="{in_file}" outdir="{out_dir}" outfile="{out_file}" exludeATWD="{exclude_atwd}" \n\n'

	with open(dag_file_name, 'a') as f:
		f.write(instructions)

	main_index+=1

in_filelist.close()