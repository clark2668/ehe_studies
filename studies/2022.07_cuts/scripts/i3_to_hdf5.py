#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/icetray-start
#METAPROJECT /home/brian/IceCube/ehe/software/build_icetray/

from icecube import icetray, dataio, dataclasses, common_variables, linefit
from I3Tray import I3Tray
from icecube import hdfwriter
import numpy as np


import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-i", type=str,
    dest="input_file",required=True,
    help="full path to the input file",
    )
parser.add_argument("-o", type=str,
    dest="output_file",required=True,
    help='''full path to the output file, without file extention. 
            That is, provide "test" not "test.i3.bz2"''',
    )
args = parser.parse_args()

filenamelist = []
filenamelist.append(args.input_file)

tray = I3Tray()
tray.AddModule("I3Reader", filenamelist=filenamelist)


tray.AddSegment(hdfwriter.I3HDFWriter, 'hdf', 
    Output=f'{args.output_file}.hdf5', 
    Keys=[
        'I3EventHeader', 'CorsikaWeightMap', 'PolyplopiaPrimary', 
        'CVMultiplicity', 'CVStatistics', 'Homogenized_QTot', 'EHELineFit',
        'LineFit'
    ],
    SubEventStreams=['InIceSplit', 'Final']
    )

tray.Execute()