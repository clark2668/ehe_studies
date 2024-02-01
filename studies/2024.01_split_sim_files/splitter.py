#!/usr/bin/env python3

import argparse, os, math

# IceCube imports
from icecube import icetray, dataio
from icecube.icetray import I3Tray

parser = argparse.ArgumentParser()
parser.add_argument("-i", type=str,
    dest="input_file", required=True,
    help="full path to the input file")
parser.add_argument("-o", type=str, 
    dest="output_dir", required=True,
    help="directory where you want things written to")
parser.add_argument("-s", type=str, 
    dest="species", required=True,
    help="species we are analyzing")

split_quantity = {
    'mu': 2,
    'nue': 8,
    'numu': 8,
    'nutau': 8,
    'tau': 10
}

def countup_qframes(in_filename):
    num_q_frames = 0
    file_in = dataio.I3File(in_filename)
    while file_in.more():
        frame = file_in.pop_frame()
        if frame.Stop == icetray.I3Frame.DAQ:
            num_q_frames+=1
    return num_q_frames

args = parser.parse_args()
input_file = args.input_file
output_dir = args.output_dir
species = args.species

if species in split_quantity:
    frames_per_file = split_quantity[species]
else:
    frames_per_file = math.floor(countup_qframes(input_file)/10.)

file_base = os.path.splitext(os.path.splitext(os.path.basename(input_file))[0])[0] # isolate the filename
print(file_base)

chunk = 0
while chunk < 10:
    
    first_q_frame = chunk * frames_per_file
    last_q_frame = first_q_frame + frames_per_file - 1
    print(f"Chunk {chunk}, FperF {frames_per_file}, First {first_q_frame}, Last {last_q_frame}")
    
    # set up output file
    out_filename = os.path.join(output_dir,
                                f'{file_base}_{chunk:01}.i3.zst'
                                )
    file_out = dataio.I3File(out_filename, 'w')
    
    # now stream over input frames
    
    file_in = dataio.I3File(input_file)
    which_q = 0
    while file_in.more():
        frame = file_in.pop_frame()
        if frame.Stop == icetray.I3Frame.DAQ:
            if which_q >= first_q_frame and which_q <= last_q_frame:
                # print(f"     Keep {which_q}")
                file_out.push(frame)
            which_q += 1
        else:
            file_out.push(frame)
        
    chunk += 1
