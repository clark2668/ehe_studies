#!/bin/bash

input_file_dir=$1
input_file_name=$2
output_file_dir=$3
output_file_name=$4

script_path=/home/brian/IceCube/ehe/ehe_studies/studies/2022.08_cal_errata_study

$script_path/redo_calib.py -i $input_file_dir/$input_file_name -o $output_file_dir/$output_file_name
