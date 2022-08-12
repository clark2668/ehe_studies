#!/bin/bash

input_file_dir=$1
input_file_name=$2
output_file_dir=$3
output_file_name=$4

script_path=/home/brian/IceCube/ehe/ehe_studies/studies/2022.07_cuts/scripts

# printf "Host is: "; hostname
# echo $_CONDOR_SCRATCH_DIR

# $script_path/i3_to_hdf5.py -i $input_file_dir/$input_file_name -o $TMPDIR/$output_file_name
$script_path/i3_to_hdf5.py -i $input_file_dir/$input_file_name -o $output_file_dir/$output_file_name
