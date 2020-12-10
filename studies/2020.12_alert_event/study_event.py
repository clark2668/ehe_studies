from icecube import icetray, dataio, dataclasses, phys_services
from I3Tray import I3Tray

from filter import HeseFilter

from icecube.phys_services.which_split import which_split
from icecube.frame_object_diff.segments import uncompress
from filter import HeseFilter, PutPulsesFromPulseMaskIntoFrame, LoopPulses

tray = I3Tray()
tray.AddModule("I3Reader", filename='alert_event.i3')

# rehydrate the I3Calibration
tray.Add(uncompress, "GCD_uncompress",
	keep_compressed=False) # will pull in baseline_gcd_134137.i3 by default

# based on https://code.icecube.wisc.edu/projects/icecube/browser/IceCube/projects/filterscripts/trunk/python/all_filters.py#L157

# first, get an appropriately named pulse series mask into the P-frame
# tray.AddModule(PutPulsesFromPulseMaskIntoFrame, 'PulseMaskPutter',
# 	pulse_mask_name='InIceDSTPulses',
# 	final_name='InIceSplitPulses',
# 	If = (which_split(split_name = 'InIceSplit')), # "check if the header.sub_event_stream == InIceSplit",
# 	Streams=[icetray.I3Frame.Physics]
# 	)

# then run the HESE filter on it
# tray.AddSegment(HeseFilter, "HeseFilter",
# 	pulses='SplitUncleanedInIcePulses',
# 	If = (which_split(split_name = 'InIceSplit')) # "check if the header.sub_event_stream == InIceSplit"
# 	)

tray.AddModule(LoopPulses, "LoopPulses",
	pulses='SplitUncleanedInIcePulses',
	Streams=[icetray.I3Frame.Physics]
	)


# tray.AddModule(TesterModule, 'tester', Streams=[icetray.I3Frame.DAQ])
# tray.Add("Dump")
tray.Add("I3Writer", filename="quick.i3.zst")
tray.Execute(5)