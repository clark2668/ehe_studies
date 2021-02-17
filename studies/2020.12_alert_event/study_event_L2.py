from icecube import icetray, dataio, dataclasses, phys_services
from I3Tray import I3Tray

from icecube.phys_services.which_split import which_split
from icecube.frame_object_diff.segments import uncompress
from filter import HeseFilter, LoopHESEPulses, LoopEHEPulses, Compare_HESE_EHE

# this starts from L2, and specifically the event file I made
# that already has the GCD file prepended

use_fadc=False
hese_pulses = 'SplitInIcePulses'

tray = I3Tray()
# tray.AddModule("I3Reader", filename='134777_8912764_gcd.i3.zst')
# tray.AddModule("I3Reader", filename='134777_8912764_L2_noFADC.i3.zst')
tray.AddModule("I3Reader", filename='134777_8912764_L2_FADC_{}.i3.zst'.format(use_fadc))

# then run the HESE filter on it
# tray.AddSegment(HeseFilter, "HeseFilter",
# 	pulses='SplitUncleanedInIcePulses',
# 	If = (which_split(split_name = 'InIceSplit')) # "check if the header.sub_event_stream == InIceSplit"
# 	)

tray.AddModule(LoopHESEPulses, "LoopPulses",
	pulses=hese_pulses,
	Streams=[icetray.I3Frame.Physics]
	)

tray.AddModule(LoopEHEPulses, "LoopPortiaPulses",
	excludeHighQE=True,
	excludeFADC=not use_fadc,
	Streams=[icetray.I3Frame.Physics]
	)

# tray.AddModule(Compare_HESE_EHE, "Compare",
# 	hese_pulses=hese_pulses,
# 	exclude_fadc= not use_fadc,
# 	table_name='comparison_overlap_{}_FADC_{}.hdf5'.format(hese_pulses, use_fadc),
# 	Streams=[icetray.I3Frame.Physics]
# 	)

# tray.Add("I3Writer", filename="quick.i3.zst")
tray.Execute()