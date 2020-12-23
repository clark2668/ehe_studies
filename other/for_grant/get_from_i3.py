# what I'm going to try and demonstrate in this code is how to loop over an I3 file
# get out the interesting information (in this case, the npe and zeniths)
# and write that to an hdf5 file that you can analyze later
# this is only one way to tackle this--we can discuss others once you have the idea of it

# python imports
import argparse
import h5py
import numpy as np

# IceCube imports
from icecube import icetray, dataio, portia, ophelia

parser = argparse.ArgumentParser()
parser.add_argument("-i", type=str, nargs='+',
	dest="input_files",
	help="full path to the input file")
parser.add_argument("-o", type=str, 
	dest="ouput_dir",
	help="directory where the output should be written to")
parser.add_argument("-v", type=int, default=0,
	dest="verbose_mode",
	help="verbosity setting; 0 = not verbose (default), 1 = verbose")

args = parser.parse_args()
input_files = args.input_files
ouput_dir = args.ouput_dir
verbose_mode = args.verbose_mode

# these are arrays which will hold our output information
npes=[]
zeniths=[]
runs=[]
event_ids=[]
subevent_ids=[]

for file in input_files:
	
	# this is one way to open the file
	file_in = dataio.I3File(file)

	i = 0
	frameId=0
	maxEvents=5e6 # big number
	while file_in.more() and i<maxEvents:
		frameId+=1
		try:
			# see if the frame we have encountered is a P-frame
			frame = file_in.pop_physics()
		except:
			continue

		# get the header for this frame
		header = frame.Get("I3EventHeader")
		run = header.run_id
		
		# skip if it's not an InIceSplit P-frame, which means it contains
		# data from the in-ice part of IceCube, not the ice-top part
		if header.sub_event_stream != "InIceSplit":
			continue 


		# so, we want to plot the number of photoelectroncs (NPE) recorded
		# and the zenith of the reconstructed particle (cos(theta))
		# the NPE is stored in the P-frame in an object called 'EHEPortiaEventSummarySRT'
		# the theta is stored in the P-frme in an object called 'EHEOpheliaParticleSRT_ImpLF'
		# so we should check that those objects are in the frame, and quit if they aren't

		# For your knowledge, "Portia" stands for PORTable Impulse Analyzer
		# and it is the tool used in the EHE analysis for figuring out how many
		# DOMS recorded signal, and also how much signal each recorded
		# it also bears the name of a Shakespeare character, which is typical
		# of the software that was designed by the Chiba group :)

		# Ophelia is the tool used in the EHE analysis for figuring out where 
		# a particle came from, so things like the particles theta/phi or azimuth/zenith
		# again, a Shakespeare naem

		if not frame.Has('EHEOpheliaParticleSRT_ImpLF'):
			continue
		if not frame.Has('EHEPortiaEventSummarySRT'):
			continue

		event_id = header.event_id
		subevent_id = header.sub_event_id

		# now, we can get the zenith from the EHEOpheliaParticleSRT_ImpLF
		zenith = frame.Get("EHEOpheliaParticleSRT_ImpLF").dir.zenith

		# and, we can get the npe from the EHEPortiaEventSummarySRT
		npe = frame.Get("EHEPortiaEventSummarySRT").GetTotalBestNPE()

		if (verbose_mode):
			print('Particle {:5}: NPE={:10.2f}, Zenith = {:.3f}'
				.format(i, npe, zenith))

		npes.append(npe)
		zeniths.append(zenith)
		runs.append(run)
		event_ids.append(event_id)
		subevent_ids.append(subevent_id)

		i+=1

# turn them all into np arrays before storing for output
npes = np.asarray(npes)
zeniths = np.asarray(zeniths)
runs = np.asarray(runs)
event_ids = np.asarray(event_ids)
subevent_ids = np.asarray(subevent_ids)


# and now, we will right them to an output file
output_file_path = "{}/output.hdf5".format(ouput_dir)
file_out = h5py.File(output_file_path, "w")

data = file_out.create_group("data")
data.create_dataset("npe",data=npes)
data.create_dataset("zenith",data=zeniths)
data.create_dataset("run",data=runs)
data.create_dataset("event_id",data=event_ids)
data.create_dataset("subevent_id",data=subevent_ids)

file_out.close()

