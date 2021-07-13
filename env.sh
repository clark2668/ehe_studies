#!/bin/bash

hostname=$HOSTNAME

if [[ $hostname == *"wisc"* ]]; then
	echo 'At UW/WIPAC filesystem! (' $hostname ')'
	export PYTHONPATH=/home/brianclark/IceCube/ehe_studies/utils/:$PYTHONPATH
	export PYTHONPATH=/home/brianclark/IceCube/modified_icetray/pip_installed_tools/:$PYTHONPATH
	export PYTHONPATH=/home/brianclark/IceCube/modified_icetray/simweights/:$PYTHONPATH
fi
#module purge
#eval `/cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/setup.sh`
#$SROOT/metaprojects/combo/V01-01-00/env-shell.sh
#/cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/RHEL_7_x86_64/metaprojects/combo/V01-00-02/env-shell.sh
#export PYTHONPATH=$PYTHONPATH:$ROOTSYS/lib
