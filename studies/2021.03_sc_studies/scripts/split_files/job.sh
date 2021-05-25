#!/bin/sh

year=$1
candle=$2
filter=$3
output_dir=$4
input_dir=$5
exclude_atwd=$6
exclude_fadc=$7

files=$input_dir/*.i3.*

export PYTHONPATH=/home/brian/IceCube/ehe_studies:$PYTHONPATH

/misc/home/brian/IceCube/ehe_studies/studies/2021.03_sc_studies/split_files.py -y $year -c $candle -f $filter -i $files -o $TMPDIR/ -atwd $exclude_atwd -fadc $exclude_fadc -r 'True'

mv $TMPDIR/*.i3.* $output_dir/.
mv $TMPDIR/*.hdf5 $output_dir/.
