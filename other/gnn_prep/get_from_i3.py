# what I'm going to try and demonstrate in this code is how to loop over an I3 file
# get out the interesting information (in this case, the energy)
# and write that to an hdf5 file that you can analyze later
# this is only one way to tackle this--we can discuss others once you have the idea of it

# python imports
import argparse
import h5py
import numpy as np

# IceCube imports
from icecube import icetray, dataio

parser = argparse.ArgumentParser()
parser.add_argument("-i", type=str, nargs='+',
	dest="input_files", required=True,
	help="full path to the input file")
parser.add_argument("-o", type=str, 
	dest="ouput_dir", required=True,
	help="directory where the output should be written to")

args = parser.parse_args()
input_files = args.input_files
ouput_dir = args.ouput_dir

# these are arrays which will hold our output information
energies = []

for file in input_files:

    print("on file {}".format(file))

    # this is one way to open the file
    file_in = dataio.I3File(file)
    
    i = 0
    frameId = 0
    maxEvents = 1E6 # some very large number

    while file_in.more() and i<maxEvents:
        print("on frame {}".format(frameId))
        frameId +=1
        try:
            # see if there's another Q frame available
            frame = file_in.pop_daq()
        except:
            continue
    
        # get the header for this frame
        header = frame.Get("I3EventHeader")
        
        # the simulation we are looking at is of muons
        # so we can get the "Monte Carlo Muon" (MCMuon) out of the frame
        # as well as extract the energy
        event = frame.Get("MCMuon")
        energy = event.energy
        energies.append(energy)
    
        i+=1 # advance counter

# turn them all into np arrays before storing for output
energies = np.asarray(energies)
# you will probably want to change the code a bit
# so that the output file name 
output_file_path = "{}/output.hdf5".format(ouput_dir)
file_out = h5py.File(output_file_path, "w")
data = file_out.create_group("data")
data.create_dataset("energies",data=energies)
file_out.close()
