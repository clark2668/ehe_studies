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


# make a replacement to write it out here at UMD
umd_out = uw_in.replace(
    "/data/exp/IceCube/",
    "/data/i3store/users/baclark/datawarehouse_copy/"
)
umd_out = f"file://{umd_out}" # do file staging
print(umd_out)

# where does it live in the grid?
# grid_in = f"http://icecube:skua@convey.icecube.wisc.edu{uw_in}"
grid_in = f"gsiftp://gridftp.icecube.wisc.edu{uw_in}"


tray = I3Tray()
# tray.context['I3FileStager'] = dataio.get_stagers() # apparently it loads the stager itself?
tray.Add(dataio.I3Reader, 'reader', FilenameList=[grid_in])

def filt(frame):
    
    filter_name = ['EHEFilter_12', 'EHEFilter_13', 'EHEAlertFilter_15',
                    'EHEAlertFilterHB_15', 'HighQFilter_17']
    passes = False
    for f in filter_name:
        if f in frame['QFilterMask']:
            if frame['QFilterMask'][f].condition_passed:
                passes = True
    return passes

tray.Add(filt)

def FilterHighCharge(frame, cut_value):
    if frame.Has('Homogenized_QTot'):
        charge = frame['Homogenized_QTot']
        if charge >= cut_value:
            return True
        else:
            return False
    else:
        return False

tray.AddModule(FilterHighCharge, 'charge cut',
               cut_value=(10**3.6),
               Streams=[icetray.I3Frame.Physics]
        )

tray.Add("I3OrphanQDropper")
tray.Add("I3Writer", FileName=umd_out)

print(time.ctime())
tray.Execute()
print(time.ctime())
