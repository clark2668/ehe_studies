from icecube import icetray, dataio, dataclasses, phys_services
from I3Tray import I3Tray


def find_event(frame, run, evt):
	keeper = False
	if frame.Has('I3EventHeader'):
		header = frame.Get('I3EventHeader')
		if header.run_id == run and header.event_id == evt:
			keeper = True
			# if it's a P frame, but not an InIceSplit, then drop the frame
			if frame.Stop == icetray.I3Frame.Physics:
				if header.sub_event_stream != "InIceSplit":
					keeper=False
	return keeper

tray = I3Tray()
run_no = 134777
evt_no = 8912764
tray.AddModule("I3Reader", filename='data/Level2_IC86.2020_data_Run00134777_Subrun00000000_00000032.i3.zst')

tray.AddModule(find_event, 'find_event', run = run_no, evt = evt_no,
	Streams=[icetray.I3Frame.Physics, icetray.I3Frame.DAQ])

tray.AddModule("I3Writer", "wrat",
	Streams=[icetray.I3Frame.Physics, icetray.I3Frame.DAQ], 
	filename=f'{run_no}_{evt_no}.i3')

tray.Execute()
tray.Finish()

# tray.Add("I3Writer", filename="quick.i3.zst")
# tray.Execute(5)