# python imports
import argparse

# IceCube imports
from icecube import icetray, dataio, dataclasses
from icecube.recclasses import I3PortiaEvent

# custom imports
import ehe_utils as ehe_utils

#/data/exp/IceCube/2011/filtered/level2pass2/1113/Run00118920/Level2pass2_IC86.2011_data_Run00118920_Subrun00000000_00000209.i3.zst

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, 
	dest="input_file",
	help="full path to the input file")

args = parser.parse_args()
input_file = args.input_file

file = dataio.I3File(input_file)

i = 0
maxEvents=5000
while file.more() and i<maxEvents:
	try:
		frame = file.pop_physics()
	except:
		continue

	# skip if it's not an InIceSplit
	if frame.Get("I3EventHeader").sub_event_stream != "InIceSplit":
		continue

	if ehe_utils.has_ehe_objects(frame):
		# print("Particle {} Has EHE objects".format(i))
		portia_npe, portia_chans = ehe_utils.get_portia_pulses_and_chans(frame)
		print("Particle {} has {} NPE and {} Chans".format(i,portia_npe, portia_chans))

	i+=1

