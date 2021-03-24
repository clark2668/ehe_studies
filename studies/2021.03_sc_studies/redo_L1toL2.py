from icecube import icetray
from icecube.icetray import I3Units
from I3Tray import I3Tray

from icecube.filterscripts.offlineL2.Rehydration import Rehydration
from icecube.filterscripts.offlineL2.level2_EHE_Calibration import EHECalibration
from icecube.filterscripts.offlineL2.level2_HitCleaning_EHE import HitCleaningEHE
from icecube.filterscripts.offlineL2.Globals import ehe_wg, ehe_wg_Qstream

from icecube.phys_services.which_split import which_split
# icetray.set_log_level(icetray.I3LogLevel.LOG_TRACE)

class drop_pframe(icetray.I3Module):
	def __init__(self, context):
		icetray.I3Module.__init__(self, context)
	# def Configure(self):
		# pass
	def Physics(self, frame):
		pass # drop P frames; leave everything else alone

tray = I3Tray()
scdata = '/misc/disk15/data/IceCube/RealData/86strings/standardcandle/2015/sc2/'
tray.AddModule("I3Reader", filenamelist=[scdata+'Level2_IC86.2014_data_Run00125920_0116_0_138_GCD.i3.gz', 
	scdata+'Level2_Run00125920_Subrun00000000_00000000.i3.bz2']) #GCD first
# tray.AddModule(drop_pframe, 'dropP')

# # do rehydration
# tray.AddSegment(Rehydration, 'rehydrator',
# 	dstfile=None,
# 	mc=False,
# 	doNotQify=False,
# 	pass2=False
# 	)

# # Create a SeededRT configuration object with the standard RT settings.
# # This object will be used by all the different SeededRT modules, i.e. the
# # modules use the same causial space and time conditions, but can use
# # different seed algorithms.
# from icecube.STTools.seededRT.configuration_services import I3DOMLinkSeededRTConfigurationService
# seededRTConfig = I3DOMLinkSeededRTConfigurationService(
# 	ic_ic_RTRadius              = 150.0*I3Units.m,
# 	ic_ic_RTTime                = 1000.0*I3Units.ns,
# 	treat_string_36_as_deepcore = False,
# 	useDustlayerCorrection      = False,
# 	allowSelfCoincidence        = True
# 	)

# tray.AddModule('I3SeededRTCleaning_RecoPulseMask_Module', 'North_seededrt',
# 	InputHitSeriesMapName  = 'SplitInIcePulses',
# 	OutputHitSeriesMapName = 'SRTInIcePulses',
# 	STConfigService        = seededRTConfig,
# 	SeedProcedure          = 'HLCCoreHits',
# 	NHitsThreshold         = 2,
# 	MaxNIterations         = 3,
# 	Streams                = [icetray.I3Frame.Physics],
# 	If = which_split(split_name='InIceSplit') & (lambda f: (ehe_wg(f)))
# 	)

# # EHE Calibration
# tray.AddSegment(EHECalibration, 'ehecalib',
# 	inPulses='CleanInIceRawData',
# 	outATWD='EHECalibratedATWD_Wave',
# 	outFADC='EHECalibratedFADC_Wave',
# 	If=lambda f: ehe_wg_Qstream(f)
# 	)

# # EHE Hit Cleaning
# tray.AddSegment(HitCleaningEHE, 'eheclean',
# 	inATWD='EHECalibratedATWD_Wave', inFADC = 'EHECalibratedFADC_Wave',
# 	# inATWD='CalibratedWaveforms_ATWD', inFADC = 'CalibratedWaveforms_FADC',
# 	If=which_split(split_name='InIceSplit') & (lambda f: ehe_wg(f))
# 	)

# tray.Add("Dump")
tray.Add("I3Writer", filename="test.i3.zst",
	# Streams=[icetray.I3Frame.TrayInfo, icetray.I3Frame.DAQ, icetray.I3Frame.Physics]
	)

tray.Execute(10)