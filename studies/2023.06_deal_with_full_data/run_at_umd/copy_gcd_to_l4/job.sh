#!/bin/bash

set=$1
run_no=$2
input_file=$3

top_dir=/data/user/brianclark/IceCube/EHE/data_full/level4_v2

cp $input_file $top_dir/$set/$run_no/.
