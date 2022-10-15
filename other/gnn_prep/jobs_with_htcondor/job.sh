# !/bin/bash

# first, load variables
inputfile=$1
outputfile=$2

# define the path to the script
script_path=/home/brianclark/IceCube/ehe_studies/other/gnn_prep

# run the script
$script_path/get_from_i3.py -i $inputfile -o $outputfile
