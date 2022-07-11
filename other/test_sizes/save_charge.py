#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/icetray-start
#METAPROJECT combo/V01-01-00

from icecube import icetray, dataio, dataclasses
from I3Tray import I3Tray
from icecube import VHESelfVeto, hdfwriter
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

# cut on the high Q filter
from icecube.filterscripts import filter_globals
def highQfilter(frame):
    if frame.Stop == icetray.I3Frame.Physics and frame.Has('FilterMask'):
        if frame['FilterMask'].get(filter_globals.HighQFilter).condition_passed:
            return 1
        else:
            return 0
    else:
        return 0

def reallyhighQfitler(frame):
    charge = frame.Get("Homogenized_QTot")
    if charge >= 4E3:
        return 1
    else:
        return 0

tray.AddModule(highQfilter, 'highQ',
    Streams=[icetray.I3Frame.Physics])

tray.AddModule(reallyhighQfitler, 'reallyhighQ',
    Streams=[icetray.I3Frame.Physics])

tray.AddModule("I3Writer", "write",
    filename=f'{args.output_file}.i3.zst',
    Streams=[icetray.I3Frame.TrayInfo, icetray.I3Frame.DAQ, icetray.I3Frame.Physics],
    DropOrphanStreams=[icetray.I3Frame.Calibration, icetray.I3Frame.DAQ]
    )

tray.Execute()
