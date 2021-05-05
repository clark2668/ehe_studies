#!/bin/bash

module purge
eval `/cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/setup.sh`
$SROOT/metaprojects/combo/V01-01-00/env-shell.sh
#/cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/RHEL_7_x86_64/metaprojects/combo/V01-00-02/env-shell.sh
#export PYTHONPATH=$PYTHONPATH:$ROOTSYS/lib
