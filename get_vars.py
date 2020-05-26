# python imports
import argparse

# IceCube imports
from icecube import icetray, dataio

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
maxEvents=1e5 # big number
while file.more() and i<maxEvents:
	try:
		frame = file.pop_physics()
	except:
		continue

	# skip if it's not an InIceSplit P-frame
	if frame.Get("I3EventHeader").sub_event_stream != "InIceSplit":
		continue

	# check if the frame contains the EHE L2 objects
	if ehe_utils.has_ehe_objects(frame):

		# if so, we get the npe and nchans as reconstructed by portia
		portia_npe, portia_chans = ehe_utils.get_portia_pulses_and_chans(frame)
		
		# and get the direction as reconstructed by ophelia
		ophelia_zenith = ehe_utils.get_ophelia_zenith(frame)
		
		# print something about the particle we found
		print("Particle {} has {} NPE, {} Chans, and {} Zenith"
			.format(i,portia_npe, portia_chans,ophelia_zenith))



	i+=1

