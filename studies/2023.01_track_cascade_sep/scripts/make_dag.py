import numpy as np
import os
from pathlib import Path

output_directory_head = '/disk20/users/brian/IceCube/juliet/level5'
# output_directory_head = '/data/i3store/users/baclark/juliet/level5'

datasets = [
            # 'numu_high_energy', 'numu_very_high_energy',
            # 'nutau_high_energy', 'nutau_very_high_energy',
            # 'tau_high_energy', 'tau_very_high_energy',
            # "21218"
            # 'nue_high_energy_nopatch', 
            # 'nue_high_energy_patch',
            # 'nue_high_energy_patch_usf'
            # 'nue_high_energy_hitmakerfix'
            # 'nue_high_energy_cmc_fine_highqcut'
            'nue_high_energy_cmc_fine'
    ]

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

        in_file = in_file.rstrip('\n') # remove the newline  
        in_dir, in_file = os.path.split(in_file) # split between path and filename
        in_mid_dir = os.path.split(in_dir)[1]
        in_file_noext = os.path.splitext(os.path.splitext(in_file)[0])[0] # strip twice (to get rid of .i3.zst)
        out_dir = output_directory_head + f'/{d}/'+ f'/{in_mid_dir}/' 
        out_file = f'{in_file_noext}'
        
        instructions = ""
        instructions += f'JOB job_{main_index} job.sub \n'
        instructions += f'VARS job_{main_index} indir="{in_dir}" infile="{in_file}" outdir="{out_dir}" outfile="{out_file}" \n\n'

        
        with open(dag_file_name, 'a') as f:
            f.write(instructions)
        main_index+=1
        
in_filelist.close()