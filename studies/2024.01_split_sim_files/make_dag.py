import argparse
import os
import numpy as np

juliet_species = ['mu', 'tau', 'nue', 'numu', 'nutau']
cor_species = [12001, 13001, 14001, 15001, 16001]

species = ['mu', 'tau', 'nue', 'numu', 'nutau']
species = [12001, 13001, 14001, 15001, 16001]
# species = ['tau', 'nue', 'numu', 'nutau']
juliet_top_dir = '/data/ana/Diffuse/EHE/simulation/juliet'
cor_top_dir = '/data/ana/Diffuse/EHE/simulation/background_corsika'

make_new_out_dirs = False
if make_new_out_dirs:
    for s in species:
        if s in juliet_species:
            for i in range(0, 12):
                folder = os.path.join(
                    juliet_top_dir,
                    s,
                    'very_high_energy',
                    'generation_split',
                    str(i)
                )
                if not os.path.exists(folder):
                    print(f"Making {folder}")
                    os.makedirs(folder)
        elif s in cor_species:
            sets = [
                '0000000-0000999', '0002000-0002999', '0004000-0004999',
                '0006000-0006999', '0008000-0008999', '0010000-0010999',
                '0001000-0001999', '0003000-0003999', '0005000-0005999',
                '0007000-0007999', '0009000-0009999'
            ]
            for set in sets:
                folder = os.path.join(
                    cor_top_dir+'_split',
                    str(s),
                    set
                )
                print(folder)
                if not os.path.exists(folder):
                    print(f"Making {folder}")
                    os.makedirs(folder)                    

get_file_list = False
if get_file_list:
    for s in species:
        if s in juliet_species:
            to_walk = os.path.join(
                juliet_top_dir,
                s,
                'very_high_energy',
                'generation'
            )
        elif s in cor_species:
            to_walk = os.path.join(
                cor_top_dir,
                str(s)
            )
        full_paths = []
        x = os.walk(to_walk)
        for root, dirnames, filenames in os.walk(to_walk):
            for filename in filenames:
                if filename.endswith('.i3.zst'):
                    full_path = os.path.join(
                        root,
                        filename
                    )
                    full_paths.append(full_path)
        full_paths = np.asarray(full_paths)
        full_paths = np.sort(full_paths)
        print(full_paths)
        np.savez(f'files_{s}.npz', full_paths=full_paths)
    
make_dag_files = True
if make_dag_files:
    for s in species:
        
        dag_filename = f"dag_{s}.dag"
        instructions = ""
        instructions += 'CONFIG config.dagman \n'
        
        
        infile = f'files_{s}.npz'
        files = np.load(infile)['full_paths']
        for i, f in enumerate(files):
            
            if i > 3000:
                break
            
            dirname = os.path.dirname(f)
            
            if s in juliet_species:
                out_dirname = dirname.replace('generation', 'generation_split')
                if not os.path.exists(out_dirname):
                    raise RuntimeError
            elif s in cor_species:
                out_dirname = dirname.replace('background_corsika', 'background_corsika_split')
                if not os.path.exists(out_dirname):
                    raise RuntimeError
            
            instructions = ""
            instructions += f'JOB job_{i} job.sub \n'
            instructions += f'VARS job_{i} infile="{f}" outdir="{out_dirname}" species="{s}" \n\n'
            
            with open(dag_filename, 'a') as fwrite:
                fwrite.write(instructions)
            
