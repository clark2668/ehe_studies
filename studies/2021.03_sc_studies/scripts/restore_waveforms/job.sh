#!/bin/bash

input_file_dir=$1
input_file_name=$2
output_file_dir=$3
output_file_name=$4

script_path=/home/brianclark/IceCube/ehe_studies/studies/2021.05_sim_studies

$script_path/save_charge.py -i $input_file_dir/$input_file_name -o $TMPDIR/$output_file_name

mv $TMPDIR/*.hdf5 $output_file_dir/.