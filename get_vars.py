# python imports
import argparse
import h5py

# IceCube imports
from icecube import icetray, dataio

# custom imports
import ehe_utils as ehe_utils # original ehe utilities
import ob_utils as ob_utils # off brand (ob) ehe utilities


# /data/exp/IceCube/2011/filtered/level2pass2a/1113/Run00118920/Level2pass2_IC86.2011_data_Run00118920_Subrun00000000_00000209.i3.zst
# update what data file we look at

parser = argparse.ArgumentParser()
parser.add_argument("-i", type=str, 
	dest="input_file",
	help="full path to the input file")
parser.add_argument("-o", type=str, 
	dest="ouput_dir",
	help="directory where the output should be written to")
parser.add_argument("-v", type=int, 
	dest="verbose_mode",
	help="verbosity setting; 0 = not verbose, 1 = verbose")

args = parser.parse_args()
input_file = args.input_file
ouput_dir = args.ouput_dir
verbose_mode = args.verbose_mode

file_in = dataio.I3File(input_file)

i = 0
maxEvents=5e2 # big number
while file_in.more() and i<maxEvents:
	try:
		frame = file_in.pop_physics()
	except:
		continue

	# skip if it's not an InIceSplit P-frame
	if frame.Get("I3EventHeader").sub_event_stream != "InIceSplit":
		continue

	# check if the frame contains the EHE L2 objects
	if ehe_utils.has_ehe_objects(frame):

		# first, get out all of the variables related to the "classic" analysis

		# get the npe and nchans as reconstructed by portia
		portia_npe, portia_chans = ehe_utils.get_portia_pulses_and_chans(frame)
		
		# get the direction as reconstructed by ophelia
		ophelia_zenith = ehe_utils.get_ophelia_zenith(frame)


		# second, get homogonized qtot, which is something we might use in the 
		# off-brand version of the analysis
		homogenized_qtot = ob_utils.get_homogenized_qtot(frame)

		if(verbose_mode):		
			print("Particle {} has {} NPE, {} Chans, and {} Zenith. Homogonized Qtot {}"
				.format(i,portia_npe, portia_chans,ophelia_zenith, homogenized_qtot))


	i+=1

output_file_path = ouput_dir+"/something.hdf5"
file_out = h5py.File(output_file_path, "w")
file_out.close()

