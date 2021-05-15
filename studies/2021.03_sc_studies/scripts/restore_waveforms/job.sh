#!/bin/bash

input_file_dir=$1
input_file_name=$2
output_file_dir=$3
output_file_name=$4
exclude_atwd=$5

export PYTHONPATH=/home/brian/IceCube/ehe_studies:$PYTHONPATH

script_path=/home/brianclark/IceCube/ehe_studies/studies/2021.03_sc_studies

$script_path/ restore_waveforms_v2.py -i $input_file_dir/$input_file_name -o $TMPDIR/$output_file_name --atwd $exclude_atwd

mv $TMPDIR/*.i3* $output_file_dir/.
mv $TMPDIR/*.hdf5* $output_file_dir/.