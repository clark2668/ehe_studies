# python imports
import argparse

# IceCube imports
from icecube import icetray, dataio, dataclasses

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
		print("Particle {} Has EHE objects".format(i))

	i+=1


# while file.more():
# 	try:
# 		frame = event_file.pop_physics()
# 	except:
# 		continue

# 	# use the real InIceSplit
# 	if frame["I3EventHeader"].sub_event_stream != "InIceSplit":
# 		continue

# 	if frame.Has("EHEOpheliaParticleSRT_ImpLF"):
# 		print("Has Ophelia Object!")

