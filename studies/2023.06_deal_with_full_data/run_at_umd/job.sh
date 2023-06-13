#!/bin/bash

input_file=$1

eval `/cvmfs/icecube.opensciencegrid.org/py3-v4.2.1/setup.sh`

env_path=$SROOT/metaprojects/icetray/v1.5.1/env-shell.sh
# script_path=/data/i3home/baclark/IceCube/ehe/ehe_studies/studies/2023.06_deal_with_full_data/run_at_umd
script_path=/home/brianclark/IceCube/ehe_studies/studies/2023.06_deal_with_full_data/run_at_umd

$env_path $script_path/highq_multi.py -i $input_file
