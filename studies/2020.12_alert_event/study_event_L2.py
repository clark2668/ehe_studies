from icecube import icetray, dataio, dataclasses, phys_services
from I3Tray import I3Tray

from icecube.phys_services.which_split import which_split
from icecube.frame_object_diff.segments import uncompress
from filter import HeseFilter, LoopHESEPulses, LoopPortiaPulses

# this starts from L2, and specifically the event file I made
# that already has the GCD file prepended

tray = I3Tray()
# tray.AddModule("I3Reader", filename='134777_8912764_gcd.i3.zst')
tray.AddModule("I3Reader", filename='134777_8912764_L2.i3.zst')

# then run the HESE filter on it
# tray.AddSegment(HeseFilter, "HeseFilter",
# 	pulses='SplitUncleanedInIcePulses',
# 	If = (which_split(split_name = 'InIceSplit')) # "check if the header.sub_event_stream == InIceSplit"
# 	)

# tray.AddModule(LoopHESEPulses, "LoopPulses",
# 	pulses='SplitInIceDSTPulses',
# 	Streams=[icetray.I3Frame.Physics]
# 	)

tray.AddModule(LoopPortiaPulses, "LoopPortiaPulses",
	Streams=[icetray.I3Frame.Physics]
	)

# tray.AddModule(TesterModule, 'tester', Streams=[icetray.I3Frame.DAQ])
# tray.Add("Dump")
tray.Add("I3Writer", filename="quick.i3.zst")
tray.Execute()