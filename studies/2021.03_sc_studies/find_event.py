from icecube import icetray, dataio, dataclasses, phys_services
from I3Tray import I3Tray


def find_event(frame, run, evts):
	keeper = False
	if frame.Has('I3EventHeader'):
		header = frame.Get('I3EventHeader')
		if header.run_id == run and header.event_id in evts:
			keeper = True
			# if it's a P frame, but not an InIceSplit, then drop the frame
			if frame.Stop == icetray.I3Frame.Physics:
				if header.sub_event_stream != "InIceSplit":
					keeper=False
	return keeper


import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-i", type=str,
	dest="input_files",required=True,
	help="full path to the input file",
	)
args = parser.parse_args()
filename = args.input_files

scdata = '/misc/disk15/data/IceCube/RealData/86strings/standardcandle/2015/sc2/'
gcdfile = 'Level2_IC86.2014_data_Run00125920_0116_0_138_GCD.i3.gz'

tray = I3Tray()
run_no = 125920
evt_no = 40325707
# evt_nos = [40088711, 40325707, 41080028, 41460269, 41662104, 42546428, 42959703, 43457772, 44548005, 44909663]
tray.AddModule("I3Reader",
	filenamelist=[filename],
	# filenamelist=[scdata+gcdfile,filename]
	)

tray.AddModule(find_event, 'find_event', run = run_no, evts = [evt_no],
	Streams=[icetray.I3Frame.Physics, icetray.I3Frame.DAQ])

tray.AddModule("I3Writer", "write",
	Streams=[icetray.I3Frame.Physics, icetray.I3Frame.DAQ], 
	filename=f'{run_no}_{evt_no}.i3.bz2')

tray.Execute()
tray.Finish()

# tray.Add("I3Writer", filename="quick.i3.zst")
# tray.Execute(5)