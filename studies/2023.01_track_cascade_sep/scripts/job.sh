#!/bin/bash

input_file_dir=$1
input_file_name=$2
output_file_dir=$3
output_file_name=$4
cor_sel=$5

eval `/cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/setup.sh` 

# umd
INSTALL_DIR='/data/i3home/baclark/IceCube/ehe/ehe_software/venv_ehe'
source $INSTALL_DIR/bin/activate
env_path=/data/i3home/baclark/IceCube/ehe/ehe_software/deps/bld_icetray/env-shell.sh
script_path=/data/i3home/baclark/IceCube/ehe/ehe_studies/studies/2023.01_track_cascade_sep/scripts


# grappa
# INSTALL_DIR='/home/brian/IceCube/ehe/ehe_software/venv_ehe'
# source $INSTALL_DIR/bin/activate
# env_path=/home/brian/IceCube/ehe/ehe_software/deps/bld_icetray/env-shell.sh
# script_path=/home/brian/IceCube/ehe/ehe_studies/studies/2023.01_track_cascade_sep/scripts

$env_path $script_path/i3_to_hdf5.py -i $input_file_dir/$input_file_name -o $output_file_dir/$output_file_name

#$script_path/i3_to_hdf5.py -i $input_file_dir/$input_file_name -o $output_file_dir/$output_file_name
