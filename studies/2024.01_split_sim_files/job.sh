#!/bin/bash

input_file=$1
output_dir=$2
species=$3

eval `/cvmfs/icecube.opensciencegrid.org/py3-v4.2.1/setup.sh`
env_path=$SROOT/metaprojects/icetray/v1.5.1/env-shell.sh
script_path=/home/brianclark/IceCube/ehe_studies/studies/2024.01_split_sim_files

$env_path $script_path/splitter.py -i $input_file  -o $output_dir -s $species
