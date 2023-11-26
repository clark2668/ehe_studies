from icecube import icetray, dataio, dataclasses, phys_services
from icecube.icetray import I3Tray

def find_event(frame, run, evt):
    keeper = False
    if frame.Has('I3EventHeader'):
        header = frame.Get('I3EventHeader')
        if header.run_id == run and header.event_id == evt:
            keeper = True
            print("Got it!")
            # # if it's a P frame, but not an InIceSplit, then drop the frame
            # if frame.Stop == icetray.I3Frame.Physics:
            # 	if header.sub_event_stream != "InIceSplit":
            # 		keeper=False
    return keeper

tray = I3Tray()
run_no = 134777
evt_no = 8912764

# # Yang's event # 2 (high signalness)
# run_no = 120117
# evt_no = 25194613
# tray.AddModule("I3Reader", filename='/data/exp/IceCube/2012/filtered/level2pass2a/0511/Run00120117/Level2pass2_IC86.2011_data_Run00120117_Subrun00000000_00000063.i3.zst')

# IC190331A (the high-E alert)
# run_no = 132379
# evt_no = 15947448
# tray.AddModule("I3Reader", filename='/data/exp/IceCube/2019/filtered/level2/0331/Run00132379/Level2_IC86.2018_data_Run00132379_Subrun00000000_00000050.i3.zst')

tray.AddModule(find_event, 'find_event', run = run_no, evt = evt_no,
    Streams=[icetray.I3Frame.Physics, icetray.I3Frame.DAQ])

tray.AddModule("I3Writer", "wrat",
    Streams=[icetray.I3Frame.Physics, icetray.I3Frame.DAQ], 
    filename=f'{run_no}_{evt_no}.i3.zst')

tray.Execute()
tray.Finish()
