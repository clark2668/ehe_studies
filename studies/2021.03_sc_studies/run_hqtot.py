#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/icetray-start
#METAPROJECT combo/V01-00-02

import argparse

from icecube import icetray
from icecube.icetray import I3Units
from I3Tray import I3Tray
from icecube import VHESelfVeto, hdfwriter

parser = argparse.ArgumentParser()
parser.add_argument("-i", type=str,
	dest="input_files",required=True,
	help="full path to the input file",
	)
parser.add_argument("-o", type=str, 
	dest="output_dir",required=True,
	help="directory where the output should be written to"
	)
args = parser.parse_args()
filename = args.input_files
output_location = args.output_dir

pulses='InIcePulses'

tray = I3Tray()
scdata = '/misc/disk15/data/IceCube/RealData/86strings/standardcandle/2015/sc2/'
tray.AddModule("I3Reader", 
	filenamelist=[scdata+'Level2_IC86.2014_data_Run00125920_0116_0_138_GCD.i3.gz', 
	scdata+filename]
	)


tray.AddModule('HomogenizedQTot', 'qtot_total', 
	Pulses=pulses,
	Output='HomogenizedQTot'
	)

tray.AddSegment(hdfwriter.I3HDFWriter, 'hdf', Output=filename+'.hdf5', 
	Keys=['I3EventHeader', 'HomogenizedQTot'], 
	SubEventStreams=['InIceSplit']
	)


tray.Add("I3Writer", filename=output_location+"/"+filename,
	Streams=[icetray.I3Frame.TrayInfo, icetray.I3Frame.DAQ, icetray.I3Frame.Physics]
	)

tray.Execute()
