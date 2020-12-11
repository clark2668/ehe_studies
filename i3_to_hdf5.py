#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/icetray-start
#METAPROJECT combo/V01-00-02


# python imports
import argparse
import h5py
import numpy as np

# IceCube imports
from icecube import icetray, dataio

# custom imports
import utils.utils_ehe as ehe_utils # original ehe utilities
import utils.utils_ob as ob_utils # off brand (ob) ehe utilities
import utils.utils_io as fh

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

npe=[]
chans=[]
zeniths=[]
homog_qtots=[]
fit_quals=[]
runs=[]
subruns=[]
parts=[]
event_ids=[]
subevent_ids=[]

for file in input_files:
	run, subrun, part = fh.get_run_subrun_part(file)
	
	print("On run %d"%run)
	file_in = dataio.I3File(file)

	i = 0
	frameId=0
	maxEvents=5e6 # big number
	while file_in.more() and i<maxEvents:
		frameId+=1
		try:
			frame = file_in.pop_physics()
		except:
			continue

		header = frame.Get("I3EventHeader") # get the header for this frame
		
		if header.sub_event_stream != "InIceSplit":
			continue # skip if it's not an InIceSplit P-frame

		if ehe_utils.has_ehe_objects(frame):

			print(frameId)
			event_id = header.event_id
			print("event id{}".format(event_id))
			subevent_id = header.sub_event_id
			portia_npe, portia_nchan = ehe_utils.get_portia_pulses_and_chans(frame)
			ophelia_zenith = ehe_utils.get_ophelia_zenith(frame)
			ophelia_fitqual = ehe_utils.get_ophelia_fitqual(frame)
			homogenized_qtot = ob_utils.get_homogenized_qtot(frame)

			if(verbose_mode):		
				print("Particle {:5}: NPE = {:10.2f} , NChan = {:3}, Zenith = {:.3f}, FitQual = {:.3f}. Homogonized Qtot = {:10.2f}"
					.format(i, portia_npe, int(portia_nchan), ophelia_zenith, ophelia_fitqual, homogenized_qtot))

			npe.append(portia_npe)
			chans.append(portia_nchan)
			zeniths.append(ophelia_zenith)
			fit_quals.append(ophelia_fitqual)
			homog_qtots.append(homogenized_qtot)

			runs.append(run)
			subruns.append(subrun)
			parts.append(part)
			event_ids.append(event_id)
			subevent_ids.append(subevent_id)

		i+=1

# turn them all into np arrays before storing for output
npe = np.asarray(npe)
chans = np.asarray(chans)
zeniths = np.asarray(zeniths)
homog_qtots = np.asarray(homog_qtots)
fit_quals = np.asarray(fit_quals)
runs = np.asarray(runs)
subruns = np.asarray(subruns)
parts = np.asarray(parts)
event_ids = np.asarray(event_ids)
subevent_ids = np.asarray(subevent_ids)


output_file_path = "{}/run{}.hdf5".format(ouput_dir, run, subrun, part)
file_out = h5py.File(output_file_path, "w")

data = file_out.create_group("data")
data.create_dataset("portia_npe",data=npe)
data.create_dataset("portia_nchan",data=chans)
data.create_dataset("ophelia_zenith",data=zeniths)
data.create_dataset("ophelia_fitqual",data=fit_quals)
data.create_dataset("homogenized_qtot",data=homog_qtots)
data.create_dataset("runs",data=runs)
data.create_dataset("subruns",data=subruns)
data.create_dataset("parts",data=parts)
data.create_dataset("event_id",data=event_ids)
data.create_dataset("subevent_id",data=subevent_ids)

file_out.close()

