from icecube import icetray, dataio, dataclasses, phys_services
from I3Tray import I3Tray
import numpy as np

def find_event_2(frame):
    keeper = False

    if frame.Has('LineFit'):
        lfvalues = frame.Get("LineFit")
        speed = lfvalues.speed
        # if (speed < 0.25) and (speed > 0.2):
        if speed < 0.15:
            keeper = True
    
    return keeper

# cut on the high Q filter
from icecube.filterscripts import filter_globals
def highQfilter(frame):
    if frame.Stop == icetray.I3Frame.Physics and frame.Has('FilterMask'):
        if frame['FilterMask'].get(filter_globals.HighQFilter).condition_passed:
            return 1
        else:
            return 0
    else:
        return 0

def reallyhighQfilter(frame):                                                                                                        
    if frame.Stop == icetray.I3Frame.Physics and frame.Has('Homogenized_QTot'):                                                      
        hqtot = frame.Get('Homogenized_QTot').value                                                                                  
        if np.log10(hqtot) < np.log10(25000):
            return 0                                                                                                                 
        else:                                                                                     
            return 1                                                                                                                 

tray = I3Tray()

filenamelist = []
# filenamelist.append('/cvmfs/icecube.opensciencegrid.org/data/GCD/GeoCalibDetectorStatus_2020.Run134142.Pass2_V0.i3.gz')
# filenamelist.append('/disk19/users/mmeier/simulations/table_based_sim/juliet/mu/high_energy/L2/1/Level2_00000999.i3.zst')

# filedir = "/data/sim/IceCube/2016/filtered/level2/neutrino-generator/21220/0000000-0000999/"
filedir = "/disk19/users/mmeier/simulations/table_based_sim/juliet/mu/high_energy/L2/1/"
from glob import glob
# files = sorted(glob(filedir + "Level2_IC86.2016_NuMu.021220.0000*.i3.zst"))
files = sorted(glob(filedir + "Level2_0000051*.i3.zst"))
for f in files:
    filenamelist.append(f)
print(filenamelist)

tray = I3Tray()
tray.AddModule("I3Reader", filenamelist=filenamelist)

tray.AddModule(highQfilter, 'highQ',
Streams=[icetray.I3Frame.Physics])

tray.AddModule(reallyhighQfilter, 'reallyhighQ',                                                                                     
    Streams=[icetray.I3Frame.Physics]                                                                                                
    )

tray.AddModule(find_event_2, 'find_event',
    Streams=[icetray.I3Frame.Physics])

tray.AddModule("I3Writer", "wrat",
    Streams=[icetray.I3Frame.Physics, icetray.I3Frame.DAQ], 
    DropOrphanStreams=[icetray.I3Frame.Calibration, icetray.I3Frame.DAQ],
    filename=f'juliet_mu_verylowspeed.i3.zst')

tray.Execute()
tray.Finish()
