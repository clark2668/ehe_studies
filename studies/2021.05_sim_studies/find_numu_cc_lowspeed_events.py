from icecube import icetray, dataio, dataclasses, phys_services
from I3Tray import I3Tray
import numpy as np

def find_event_2(frame):
    keeper = False
    if frame.Has('I3MCWeightDict'):
        weight_dit = frame.Get('I3MCWeightDict')
        int_type = int(weight_dit['InteractionType'])
           
        if frame.Has('HitStatisticsValues'):
            hsvalues = frame.Get('HitStatisticsValues')
            cogz = hsvalues.cog.z

            if frame.Has('LineFit'):
                lfvalues = frame.Get("LineFit")
                speed = lfvalues.speed

            # if frame.Has('Homogenized_QTot'):
                # hqtot = frame.Get('Homogenized_QTot').value

                if(int_type < 2 and speed < 0.2):
                    print('Int Type {}, COG z {}, Speed {}'.format(int_type, cogz, speed))
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
        if np.log10(hqtot) < 4.4:                                                                                                      
            return 0                                                                                                                 
        else:                                                                                                                        
            return 1                                                                                                                 

tray = I3Tray()


# 

filenamelist = []
filenamelist.append('/cvmfs/icecube.opensciencegrid.org/data/GCD/GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE_withScaledNoise.i3.gz')
filenamelist.append('/data/sim/IceCube/2016/filtered/level2/neutrino-generator/21220/0000000-0000999/Level2_IC86.2016_NuMu.021220.000000.i3.zst')
filenamelist.append('/data/sim/IceCube/2016/filtered/level2/neutrino-generator/21220/0000000-0000999/Level2_IC86.2016_NuMu.021220.000001.i3.zst')
filenamelist.append('/data/sim/IceCube/2016/filtered/level2/neutrino-generator/21220/0000000-0000999/Level2_IC86.2016_NuMu.021220.000002.i3.zst')
filenamelist.append('/data/sim/IceCube/2016/filtered/level2/neutrino-generator/21220/0000000-0000999/Level2_IC86.2016_NuMu.021220.000003.i3.zst')
filenamelist.append('/data/sim/IceCube/2016/filtered/level2/neutrino-generator/21220/0000000-0000999/Level2_IC86.2016_NuMu.021220.000004.i3.zst')
filenamelist.append('/data/sim/IceCube/2016/filtered/level2/neutrino-generator/21220/0000000-0000999/Level2_IC86.2016_NuMu.021220.000005.i3.zst')

# filedir = "/data/sim/IceCube/2016/filtered/level2/neutrino-generator/21220/0000000-0000999/"
# from glob import glob
# files = sorted(glob(filedir + "Level2_IC86.2016_NuMu.021220.0000*.i3.zst"))
# for f in files:
#     filenamelist.append(f)

tray = I3Tray()
tray.AddModule("I3Reader", filenamelist=filenamelist)

tray.AddModule(highQfilter, 'highQ',
Streams=[icetray.I3Frame.Physics])

tray.AddModule(reallyhighQfilter, 'reallyhighQ',                                                                                     
    Streams=[icetray.I3Frame.Physics]                                                                                                
    )

# and also need things like the COG (hit statistics)
from icecube.common_variables import hit_statistics
hsoutput_name = 'HitStatisticsValues'
pulses = 'SplitInIcePulses'
tray.AddSegment(hit_statistics.I3HitStatisticsCalculatorSegment, 'hs',
    PulseSeriesMapName=pulses,
    OutputI3HitStatisticsValuesName=hsoutput_name
)
tray.AddModule(find_event_2, 'find_event',
    Streams=[icetray.I3Frame.Physics])

tray.AddModule("I3Writer", "wrat",
    Streams=[icetray.I3Frame.Physics, icetray.I3Frame.DAQ], 
    DropOrphanStreams=[icetray.I3Frame.Calibration, icetray.I3Frame.DAQ],
    filename=f'numu_cc_lowspeed.i3.zst')

tray.Execute()
tray.Finish()

# tray.Add("I3Writer", filename="quick.i3.zst")
# tray.Execute(5)