#!/bin/sh

year=$1
candle=$2
filter=$3
output_dir=$4

files=/disk19/users/brian/IceCube/standard_candle/2015/w_hqtot/Level2_Run00125920_Subrun00000000_00000*.i3.bz2

/misc/home/brian/IceCube/ehe_studies/studies/2021.03_sc_studies/split_files.py -y $year -c $candle -f $filter -i $files -o $TMPDIR/

mv $TMPDIR/*.i3.* $output_dir/.
