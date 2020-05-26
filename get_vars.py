# python imports
import argparse
import h5py
import numpy as np

# IceCube imports
from icecube import icetray, dataio

# custom imports
import ehe_utils as ehe_utils # original ehe utilities
import ob_utils as ob_utils # off brand (ob) ehe utilities
import file_handler as fh


# /data/exp/IceCube/2011/filtered/level2pass2a/1113/Run00118920/Level2pass2_IC86.2011_data_Run00118920_Subrun00000000_00000209.i3.zst
# update what data file we look at

parser = argparse.ArgumentParser()
parser.add_argument("-i", type=str, 
	dest="input_file",
	help="full path to the input file")
parser.add_argument("-o", type=str, 
	dest="ouput_dir",
	help="directory where the output should be written to")
parser.add_argument("-v", type=int, default=0,
	dest="verbose_mode",
	help="verbosity setting; 0 = not verbose (default), 1 = verbose")

args = parser.parse_args()
input_file = args.input_file
ouput_dir = args.ouput_dir
verbose_mode = args.verbose_mode

run, subrun, part = fh.get_run_subrun_part(input_file)
file_in = dataio.I3File(input_file)

npe=[]
chans=[]
zeniths=[]
homog_qtots=[]
event_ids=[]
subevent_ids=[]

i = 0
maxEvents=5e5 # big number
while file_in.more() and i<maxEvents:
	try:
		frame = file_in.pop_physics()
	except:
		continue

	header = frame.Get("I3EventHeader") # get the header for this frame
	
	if header.sub_event_stream != "InIceSplit":
		continue # skip if it's not an InIceSplit P-frame

	# check if the frame contains the EHE L2 objects
	if ehe_utils.has_ehe_objects(frame):

		# first, get some bookkeeping information
		event_id = header.event_id
		subevent_id = header.sub_event_id

		# then, get out all of the variables related to the "classic" analysis

		# npe and nchans as reconstructed by portia
		portia_npe, portia_nchan = ehe_utils.get_portia_pulses_and_chans(frame)
		
		# direction as reconstructed by ophelia
		ophelia_zenith = ehe_utils.get_ophelia_zenith(frame)


		# then, get variables we might use in an "offbrand" analysis
		homogenized_qtot = ob_utils.get_homogenized_qtot(frame)

		if(verbose_mode):		
			print("Particle {:5}: NPE = {:10.2f} , NChan = {:3}, Zenith = {:.3f}. Homogonized Qtot = {:10.2f}"
				.format(i, portia_npe, int(portia_nchan), ophelia_zenith, homogenized_qtot))


		# store everything
		npe.append(portia_npe)
		chans.append(portia_nchan)
		zeniths.append(ophelia_zenith)
		homog_qtots.append(homogenized_qtot)
		event_ids.append(event_id)
		subevent_ids.append(subevent_id)

	i+=1

# turn them all into np arrays before storying for output
npe = np.asarray(npe)
chans = np.asarray(chans)
zeniths = np.asarray(zeniths)
homog_qtots = np.asarray(homog_qtots)
event_ids = np.asarray(event_ids)
subevent_ids = np.asarray(subevent_ids)

# write information out to hdf5 file

output_file_path = "{}/Run{}_Subrun{}_Part{}.hdf5".format(ouput_dir, run, subrun, part)
file_out = h5py.File(output_file_path, "w")

data = file_out.create_group("data")
data.create_dataset("portia_npe",data=npe)
data.create_dataset("portia_nchan",data=npe)
data.create_dataset("ophelia_zenith",data=zeniths)
data.create_dataset("homogenized_qtot",data=homog_qtots)
data.create_dataset("event_id",data=event_ids)
data.create_dataset("subevent_id",data=subevent_ids)

metadata = file_out.create_group("metadata")
metadata.attrs['run_id'] = run

file_out.close()

