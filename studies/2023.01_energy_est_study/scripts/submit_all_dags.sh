#!/bin/bash

FILES=dagman_*.dag
for f in $FILES
do
	condor_submit_dag $f
done