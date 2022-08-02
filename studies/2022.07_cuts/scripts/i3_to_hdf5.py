#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/icetray-start
#METAPROJECT combo/V01-01-00

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
tray.context['I3FileStager'] = dataio.get_stagers()
tray.AddModule("I3Reader", filenamelist=filenamelist)

# high Q filter (must have >1000 PE)
def highQfilter(frame):
    try:
        charge = frame.Get("Homogenized_QTot")
        if charge >= 1E3:
            return 1
        else:
            return 0
    except:
        return 0

tray.AddModule(highQfilter, 'highQ',
    Streams=[icetray.I3Frame.Physics],
    )

tray.AddSegment(hdfwriter.I3HDFWriter, 'hdf', 
    Output=f'{args.output_file}.hdf5', 
    Keys=[
        'I3EventHeader', 'CorsikaWeightMap', 'PolyplopiaPrimary', 'I3MCWeightDict',
        'CVMultiplicity', 'CVStatistics', 'Homogenized_QTot', 'EHELineFit',
        'LineFit'
    ],
    SubEventStreams=['InIceSplit', 'Final']
    )

tray.Execute()