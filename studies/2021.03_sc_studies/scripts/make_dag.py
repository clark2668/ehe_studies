dag_file_name = 'dagman_hqtot.dag'
instructions = ""
instructions += 'CONFIG config.dagman\n\n'

with open(dag_file_name, 'w') as f:
	f.write(instructions)

outputdir="/misc/disk19/users/brian/IceCube/standard_candle/2015"

input_file = open('list_2015_sc2.txt')
for i, line in enumerate(input_file):
	instructions = ""
	instructions += f'JOB job_{i} job.sub\n'
	instructions += f'VARS job_{i} file="{line[:-1]}" outputdir="{outputdir}" \n\n'

	with open(dag_file_name, 'a') as f:
		f.write(instructions)
