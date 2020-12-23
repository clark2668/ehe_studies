# all I'm trying to demonstrate here is how to get the data *back out* of the hdf5 file
# which you can then use to make plots

import h5py
import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("-f", type=str, nargs='+',
	dest="input_files",
	help="paths to input files (absolute or relative)")

args = parser.parse_args()
files = args.input_files

for file in files:
	print(file)
	file_in = h5py.File(file, "r")
	data = file_in['data']
	npes = np.asarray(data['npe'])
	zeniths = np.asarray(data['zenith'])
	print(npes)

	# now, you have a numpy array containing the npes and zeniths
	# and you can use them to make plots

	file_in.close()
