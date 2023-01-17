import numpy as np
import os
from pathlib import Path

# output_directory_head = '/disk20/users/brian/IceCube/corsika/level4/hdf5'
# output_directory_head = '/disk20/users/brian/IceCube/juliet/millipede_study'
# output_directory_head = '/disk20/users/brian/IceCube/nugen/hdf5/'
output_directory_head = '/disk20/users/brian/IceCube/corsika/millipede_study'


datasets = ['nue_high_energy', 'numu_high_energy', 
            'nutau_high_energy', 'mu_high_energy', 
            'tau_high_energy']
datasets = ['nue_very_high_energy', 'numu_very_high_energy', 
            'nutau_very_high_energy', 'mu_very_high_energy', 
            'tau_very_high_energy']
datasets = [21962, 22023]
datasets = ['numu_high_energy']
datasets = [22023]

for d in datasets:
    print(d)
        
    dag_file_name = f'dagman_{d}.dag'
    instructions = ""
    instructions += 'CONFIG config.dagman\n'
    instructions += "\n\n"
    with open(dag_file_name, 'w') as f:
        f.write(instructions)

    main_index = 0
    
    in_filelist = open(f'files_{d}.txt', 'r')
    for in_file in in_filelist.readlines():

        # if main_index > 2:
            # continue

        in_file = in_file.rstrip('\n') # remove the newline  
        in_dir, in_file = os.path.split(in_file) # split between path and filename
        in_mid_dir = os.path.split(in_dir)[1]
        in_file_noext = os.path.splitext(os.path.splitext(in_file)[0])[0] # strip twice (to get rid of .i3.zst)
        out_dir = output_directory_head + f'/{d}/'+ f'/{in_mid_dir}/' 
        out_file = f'{in_file_noext}'
        
        # if d not in [21962, 22023]:
        if 1==1:
            instructions = ""
            instructions += f'JOB job_{main_index} job.sub \n'
            instructions += f'VARS job_{main_index} indir="{in_dir}" infile="{in_file}" outdir="{out_dir}" outfile="{out_file}" \n\n'
        # else:
        #     # special
        #     # in the case of EHE corsika, we want to write out two files
        #     # one for proton, and one for iron
        #     instructions = ""
        #     instructions += f'JOB job_{main_index} job.sub \n'
        #     instructions += f'VARS job_{main_index} indir="{in_dir}" infile="{in_file}" outdir="{out_dir}" outfile="{out_file}_p" corsel="p"\n\n'
        #     main_index+=1
            
        #     instructions += ""
        #     instructions += f'JOB job_{main_index} job.sub \n'
        #     instructions += f'VARS job_{main_index} indir="{in_dir}" infile="{in_file}" outdir="{out_dir}" outfile="{out_file}_fe" corsel="fe"\n\n'

        
        with open(dag_file_name, 'a') as f:
            f.write(instructions)
        main_index+=1
        
in_filelist.close()