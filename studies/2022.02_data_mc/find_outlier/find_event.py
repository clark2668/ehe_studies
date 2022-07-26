from icecube import icetray, dataio, dataclasses, phys_services
from I3Tray import I3Tray

def find_event(frame, run, evt):
    keeper = False
    if frame.Has('I3EventHeader'):
        header = frame.Get('I3EventHeader')
        if header.run_id == run and header.event_id == evt:
            keeper = True
            # # if it's a P frame, but not an InIceSplit, then drop the frame
            # if frame.Stop == icetray.I3Frame.Physics:
            # 	if header.sub_event_stream != "InIceSplit":
            # 		keeper=False
    return keeper

import glob

run_no = 119510
evt_no = 19598412
file_list = glob.glob("/data/exp/IceCube/2012/filtered/level2pass2a/0130/Run00119510/*_00000050.i3.zst")

run_no = 120680
evt_no = 1605018
file_list = glob.glob("/data/exp/IceCube/2012/filtered/level2pass2a/0916/Run00120680/*_00000004.i3.zst")

run_no = 122590
evt_no = 45181927
file_list = glob.glob("/data/exp/IceCube/2013/filtered/level2pass2a/0624/Run00122590/*_00000121.i3.zst")

run_no = 123640
evt_no = 1035095
file_list = glob.glob("/data/exp/IceCube/2014/filtered/level2pass2a/0102/Run00123640/*_00000002.i3.zst")



tray = I3Tray()
# # run_no = 78700221
# # evt_no = 4140
# # tray.AddModule("I3Reader", 
# #     filename='/data/sim/IceCube/2016/filtered/level2/CORSIKA-in-ice/20787/0000000-0000999/Level2_IC86.2016_corsika.020787.000220.i3.zst'
# #     )

# # run_no = 78700509
# # evt_no = 1439
# # tray.AddModule("I3Reader", 
# #     filename='/data/sim/IceCube/2016/filtered/level2/CORSIKA-in-ice/20787/0000000-0000999/Level2_IC86.2016_corsika.020787.000508.i3.zst'
# #     )

tray.AddModule("I3Reader", filenamelist=file_list )

tray.AddModule(find_event, 'find_event', run = run_no, evt = evt_no,
    Streams=[icetray.I3Frame.Physics, icetray.I3Frame.DAQ])

tray.AddModule("I3Writer", "write",
    Streams=[icetray.I3Frame.Physics, icetray.I3Frame.DAQ], 
    filename=f'{run_no}_{evt_no}.i3.zst')

tray.Execute()
tray.Finish()
