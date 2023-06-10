#!/usr/bin/env python3

import sys, time, argparse
from glob import glob
from icecube import dataio, phys_services, icetray
from I3Tray import I3Tray

parser = argparse.ArgumentParser()
parser.add_argument("-i", type=str, 
    dest="infile",
    help="which dataset to read")

args = parser.parse_args()
infile = args.infile

uw_in = infile

remote_transfer = True
if remote_transfer:
    
    # where does it live in the grid?
    # grid_in = f"http://icecube:skua@convey.icecube.wisc.edu{uw_in}"
    the_in_file = f"gsiftp://gridftp.icecube.wisc.edu{uw_in}"
    
    # make a replacement to write it out here at UMD
    the_out_file = uw_in.replace(
        "/data/exp/IceCube/",
        "/data/i3store/users/baclark/datawarehouse_copy/"
    )
else:
    
    the_in_file = uw_in
    
    # put it in Brian's copy of the datawarehouse
    the_out_file = uw_in.replace(
        "/data/exp/IceCube/",
        "/data/user/brianclark/IceCube/EHE/datawarehouse_copy/"
    )

# the_out_file = f"file://{the_out_file}" # do file staging
print(f"The out file is: {the_out_file}")


tray = I3Tray()
# tray.context['I3FileStager'] = dataio.get_stagers() # apparently it loads the stager itself?
tray.Add(dataio.I3Reader, 'reader', FilenameList=[the_in_file])


def filt(frame):
    
    filter_name = ['EHEFilter_12', 'EHEFilter_13', 'EHEAlertFilter_15',
                    'EHEAlertFilterHB_15', 'HighQFilter_17']
    passes = False
    for f in filter_name:
        if f in frame['QFilterMask']:
            if frame['QFilterMask'][f].condition_passed:
                passes = True
    return passes


def FilterHighCharge(frame, cut_value):
    if frame.Has('Homogenized_QTot'):
        charge = frame['Homogenized_QTot']
        if charge >= cut_value:
            return True
        else:
            return False
    else:
        return False


moving_data = False
if moving_data:

    tray.Add(filt)

    tray.AddModule(FilterHighCharge, 'charge cut',
                cut_value=(10**3.6),
                Streams=[icetray.I3Frame.Physics]
            )

    tray.Add("I3OrphanQDropper")

tray.Add("I3Writer", FileName=the_out_file)

print(time.ctime())
tray.Execute()
print(time.ctime())
