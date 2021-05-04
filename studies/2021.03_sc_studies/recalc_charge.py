#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/icetray-start
#METAPROJECT /home/brian/IceCube/modified_icetray/bld/

from icecube import icetray, dataio, dataclasses
from I3Tray import I3Tray
from icecube import VHESelfVeto, hdfwriter

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-i", type=str,
	dest="input_files",required=True,
	help="full path to the input file",
	)
parser.add_argument("-o", type=str, 
	dest="output_dir",required=True,
	help="directory where the output should be written to"
	)
parser.add_argument("-y", type=int,
	dest="year", default=2015,
	help="Which year of standard candle, e.g. 2015"
	)
parser.add_argument("-c", type=int,
	dest="candle", default=2,
	help="Which candle, e.g. 2"
	)
parser.add_argument("-f", type=int,
	dest="filter_setting", required=True,
	help="Which standard candle filter setting, e.g. 0, 1, 2, ..."
	)

args = parser.parse_args()
filename = args.input_files


tray = I3Tray()
scdata = '/misc/disk15/data/IceCube/RealData/86strings/standardcandle/2015/sc2/'
tray.AddModule("I3Reader", filenamelist=[scdata+'Level2_IC86.2014_data_Run00125920_0116_0_138_GCD.i3.gz', 
	args.input_files]) #GCD first

pulse_name='SplitInIcePulses' # use SPLITInIcePulses

tray.AddModule('HomogenizedQTot', 'qtot_total',
	Pulses=pulse_name,
	Output='HomogenziedQtot_SplitInIcePulses'
	)

tray.AddSegment(hdfwriter.I3HDFWriter, 'hdf', 
	Output="{}/y{}_c{}_f{}_hqtot_splitinnicepulses.hdf5".format(args.output_dir, args.year, args.candle, args.filter_setting), 
	Keys=['I3EventHeader', 'HomogenizedQTot', 'HomogenziedQtot_SplitInIcePulses', 'EHEPortiaEventSummarySRT'], 
	SubEventStreams=['InIceSplit']
	)

tray.Add("I3Writer", 
	filename="{}/y{}_c{}_f{}_hqtot_splitinnicepulses.i3.bz2".format(args.output_dir, args.year, args.candle, args.filter_setting),
	Streams=[icetray.I3Frame.DAQ, icetray.I3Frame.Physics]
	)

tray.Execute()
