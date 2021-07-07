#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/icetray-start
#METAPROJECT combo/V01-01-00

from icecube import icetray, dataio, dataclasses
from I3Tray import I3Tray
from icecube import VHESelfVeto, hdfwriter

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
args = parser.parse_args()


tray = I3Tray()
tray.AddModule("I3Reader", filenamelist=[args.input_file]) #GCD first

def has_needed(frame):
	return frame.Has('EHEPortiaEventSummarySRT') and frame.Has('EHEOpheliaParticleSRT_ImpLF')

# drop every frame that doesn't have the information
tray.AddModule(has_needed, 'something',
	Streams=[icetray.I3Frame.Physics],
	)

tray.AddSegment(hdfwriter.I3HDFWriter, 'hdf', 
	Output=f'{args.output_file}.hdf5', 
	Keys=['I3MCWeightDict', 'I3EventHeader', 'I3MCTree_preMuonProp', 'PolyplopiaPrimary',
    'I3PrimaryInjectorInfo', 'I3CorsikaWeight',
    'Homogenized_QTot', 'LineFit', 
    'EHEPortiaEventSummarySRT', 'EHEOpheliaParticleSRT_ImpLF'], 
	SubEventStreams=['InIceSplit']
	)

if args.save_i3file:
	tray.AddModule("I3Writer", "write",
		filename=f'{args.output_file}.i3.zst',
		Streams=[icetray.I3Frame.DAQ, icetray.I3Frame.Physics],
		DropOrphanStreams=[icetray.I3Frame.Calibration, icetray.I3Frame.DAQ]
		)

tray.Execute()
