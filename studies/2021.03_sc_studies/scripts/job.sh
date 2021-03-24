#!/bin/sh

file=$1
output_dir=$2

/misc/home/brian/IceCube/ehe_studies/studies/2021.03_sc_studies/run_hqtot.py -i $file -o $TMPDIR/

mv $TMPDIR/*.i3.bz2* $output_dir/.
