#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/icetray-start
#METAPROJECT combo/V01-01-00

from icecube import icetray, dataio, dataclasses
from I3Tray import I3Tray
from icecube import VHESelfVeto, hdfwriter
import numpy as np

import sys
sys.path.append('/home/brianclark/IceCube/ehe_studies/utils/')

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
parser.add_argument("-s", type=bool,
    dest="save_i3file",required=False, default=False,
    help='should save i3 file',
    )
parser.add_argument("-g", type=str,
    dest="gcd_file",required=False, default=None,
    help='gcd file',
    )
# /cvmfs/icecube.opensciencegrid.org/data/GCD/GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE_withScaledNoise.i3.gz
args = parser.parse_args()

filenamelist = []
if args.gcd_file is not None:
    filenamelist.append(args.gcd_file)
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

tray.AddModule(highQfilter, 'highQ',
    Streams=[icetray.I3Frame.Physics])


def reallyhighQfilter(frame):
    if frame.Stop == icetray.I3Frame.Physics and frame.Has('Homogenized_QTot'):
        hqtot = frame.Get('Homogenized_QTot').value
        if np.log10(hqtot) < 4:
            return 0
        else:
            return 1
tray.AddModule(reallyhighQfilter, 'reallyhighQ',
    Streams=[icetray.I3Frame.Physics]
    )

# tray.AddModule('Dump')
# def printy(frame, thing):
#     print(frame.Get(thing))
# tray.AddModule(printy, thing='EHEOpheliaParticleSRT_ImpLF')

from icecube import linefit
tray.AddSegment(linefit.simple, 'yo', inputResponse='SRTInIcePulses',
    fitName='yo'
    )

from utils_pulses import get_linefit_quality
tray.AddModule(get_linefit_quality, 'LFqual', 
    linefit_name='LineFit',
    linefit_params_name='LineFitParams',
    # pulses_name='SplitInIcePulses',
    pulses_name='SRTInIcePulses',
    # pulses_name='yo_debiasedPulses',
    output_name='LineFitQuality',
    )




tray.AddSegment(hdfwriter.I3HDFWriter, 'hdf', 
    Output=f'{args.output_file}_LF.hdf5', 
    Keys=['I3MCWeightDict', 'I3EventHeader', 'PolyplopiaPrimary', 'CorsikaWeightMap',
    'InteractingNeutrino', 'PrimaryNeutrino', 'VertexPosition',
    'Homogenized_QTot', 'LineFit', 'CascadeFillRatio_L3',
    'EHEPortiaEventSummarySRT', 'EHEOpheliaParticleSRT_ImpLF', 'EHEOpheliaSRT_ImpLF'], 
    SubEventStreams=['InIceSplit']
    )

if args.save_i3file:
    tray.AddModule("I3Writer", "write",
        filename=f'{args.output_file}_LF.i3.zst',
        Streams=[icetray.I3Frame.TrayInfo, icetray.I3Frame.DAQ, icetray.I3Frame.Physics],
        DropOrphanStreams=[icetray.I3Frame.Calibration, icetray.I3Frame.DAQ]
        )

tray.Execute()
