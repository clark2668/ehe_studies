#!/usr/bin/env python3

import sys, time, argparse, glob, copy, os
import glob
from icecube import dataio, phys_services, icetray
from I3Tray import I3Tray

parser = argparse.ArgumentParser()
parser.add_argument("-i", type=str, 
    dest="indir",
    help="directory to read")

args = parser.parse_args()
indir = args.indir

# put it in Brian's copy of the datawarehouse
outdir = indir.replace(
    "/data/exp/IceCube/",
    # "file:///data/user/brianclark/IceCube/EHE/datawarehouse_copy/"
    "/data/user/brianclark/IceCube/EHE/datawarehouse_copy/"
    )
print(f"The out dir is: {outdir}")

file_list = sorted(glob.glob(f"{indir}/*.i3.zst"))

# get rid of the GCD file
new_file_list = []
for f in file_list:
    if 'GCD' not in f:
        # new_file_list.append(f"file://{f}")
        new_file_list.append(f)
file_list = copy.deepcopy(new_file_list)
del new_file_list # careful for cleanup
print(f" Looping over {len(file_list)} files")

# now work out the file pattern
demo_file_path = file_list[0]
demo_file = os.path.split(demo_file_path)[1] # the second entry
splitted = demo_file.split('_Subrun')[0]
out_file_format = f"{splitted}_merged_%08d.i3.gz"


tray = I3Tray()
tray.context['I3FileStager'] = dataio.get_stagers()
tray.Add("I3Reader", 'reader', FilenameList=file_list)


def filt(frame):
    filter_name = ['EHEFilter_12', 'EHEFilter_13', 'EHEAlertFilter_15',
                    'EHEAlertFilterHB_15', 'HighQFilter_17']
    passes = False
    for f in filter_name:
        if f in frame['QFilterMask'] and not passes:
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


moving_which = 'data'
if moving_which == 'data':

    tray.AddModule(filt)

    tray.AddModule(FilterHighCharge, 'charge cut',
                cut_value=(10**4),
                Streams=[icetray.I3Frame.Physics]
            )

    tray.Add("I3OrphanQDropper")
    
    keeper_streams = [
                    #   icetray.I3Frame.TrayInfo, 
                      icetray.I3Frame.DAQ,
                      icetray.I3Frame.Physics]

elif moving_which == 'gcd':
    
    keeper_streams = [
                    # icetray.I3Frame.TrayInfo, 
                    icetray.I3Frame.Geometry,
                    icetray.I3Frame.Calibration, icetray.I3Frame.DetectorStatus]


tray.Add("I3MultiWriter",
         Filename=f"{outdir}/{out_file_format}",
         SizeLimit=10**6,
         Streams=keeper_streams,
        #  MetadataStreams=[icetray.I3Frame.TrayInfo]
         )

print(time.ctime())
tray.Execute()
print(time.ctime())

usagemap = tray.Usage()
for mod in usagemap:
    print(mod)