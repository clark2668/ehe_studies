#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/icetray-start
#METAPROJECT /home/brian/IceCube/modified_icetray/bld/

import argparse

from icecube import icetray, dataio, dataclasses, hdfwriter, phys_services
from I3Tray import I3Tray
import tools

parser = argparse.ArgumentParser()
parser.add_argument("-i", type=str, nargs='+',
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
output_location = args.output_dir
year = args.year
candle = args.candle
filter = args.filter_setting

pulses='InIcePulses'

tray = I3Tray()

filenamelist = []
for file in filename:
	filenamelist.append(file)

tray.AddModule("I3Reader",  
	filenamelist=filenamelist 
	)

start, stop, qmin, qmax = tools.get_start_stop(year, candle, filter)
qmin = -10
qmax = 1E10
# make the qmin and qmax so that every event possible will pass the cut
start = dataclasses.I3Time(start)
stop = dataclasses.I3Time(stop)

tray.AddModule(tools.cut_by_config, 'cut',
	start=start,stop=stop,qmin=qmin,qmax=qmax,
	Streams=[icetray.I3Frame.Physics]
	)

tray.Add("I3OrphanQDropper")

tray.AddSegment(hdfwriter.I3HDFWriter, 'hdf', 
	Output=output_location+"/"+'y{}_c{}_f{}.hdf5'.format(year,candle,filter), 
	Keys=['I3EventHeader', 'HomogenizedQTot', 'HomogenizedQTot_DeepMagSix',
	'PortiaEventSummarySRT', 'PortiaEventSummarySRT_DeepMagSix'], 
	SubEventStreams=['InIceSplit'],
	)

tray.Add("I3Writer", 
	filename=output_location+"/"+'y{}_c{}_f{}.i3.zst'.format(year,candle,filter),
	Streams=[icetray.I3Frame.DAQ, icetray.I3Frame.Physics]
	)

tray.Execute()
