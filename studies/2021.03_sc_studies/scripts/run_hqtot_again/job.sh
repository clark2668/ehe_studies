#!/bin/sh

year=$1
candle=$2
filter=$3
output_dir=$4

files="/misc/disk19/users/brian/IceCube/standard_candle/2015/step3_addwave/y${year}_c${candle}_f${filter}_waves.i3.bz2"

/misc/home/brian/IceCube/ehe_studies/studies/2021.03_sc_studies/recalc_charge.py -y $year -c $candle -f $filter -i $files -o $TMPDIR/

mv $TMPDIR/*.i3.* $output_dir/.
mv $TMPDIR/*hdf5 $output_dir/.

