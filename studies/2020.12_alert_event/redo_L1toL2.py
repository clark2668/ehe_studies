from icecube import icetray, dataio, dataclasses, phys_services
from icecube.icetray import I3Units
from I3Tray import I3Tray

from icecube.filterscripts.offlineL2.Rehydration import Rehydration
from icecube.filterscripts.offlineL2.level2_EHE_Calibration import EHECalibration
from icecube.filterscripts.offlineL2.level2_HitCleaning_EHE import HitCleaningEHE
from icecube.filterscripts.offlineL2.Globals import ehe_wg, ehe_wg_Qstream

from icecube.phys_services.which_split import which_split

tray = I3Tray()
tray.AddModule("I3Reader", filenamelist=['data/Level2_IC86.2020_data_Run00134777_1209_80_564_GCD.i3.zst', '134777_8912764_L1.i3']) #GCD first

# we are following closely the official L2 scripts
# https://github.com/icecube/IceTrayCombo/blob/d4cff8945bb5f369bc4268b0c6189881129d235c/filterscripts/python/offlineL2/level2_all_filters.py

# do rehydration
tray.AddSegment(Rehydration, 'rehydrator',
	dstfile=None,
	mc=False,
	doNotQify=False,
	pass2=False
	)

# Create a SeededRT configuration object with the standard RT settings.
# This object will be used by all the different SeededRT modules, i.e. the
# modules use the same causial space and time conditions, but can use
# different seed algorithms.
from icecube.STTools.seededRT.configuration_services import I3DOMLinkSeededRTConfigurationService
seededRTConfig = I3DOMLinkSeededRTConfigurationService(
	ic_ic_RTRadius              = 150.0*I3Units.m,
	ic_ic_RTTime                = 1000.0*I3Units.ns,
	treat_string_36_as_deepcore = False,
	useDustlayerCorrection      = False,
	allowSelfCoincidence        = True
	)

tray.AddModule('I3SeededRTCleaning_RecoPulseMask_Module', 'North_seededrt',
	InputHitSeriesMapName  = 'SplitInIcePulses',
	OutputHitSeriesMapName = 'SRTInIcePulses',
	STConfigService        = seededRTConfig,
	SeedProcedure          = 'HLCCoreHits',
	NHitsThreshold         = 2,
	MaxNIterations         = 3,
	Streams                = [icetray.I3Frame.Physics],
	If = which_split(split_name='InIceSplit') & (lambda f: (ehe_wg(f)))
	)

# EHE Calibration
tray.AddSegment(EHECalibration, 'ehecalib',
	inPulses='CleanInIceRawData',
	outATWD='EHECalibratedATWD_Wave',
	outFADC='EHECalibratedFADC_Wave',
	If=lambda f: ehe_wg_Qstream(f)
	)

# EHE Hit Cleaning
tray.AddSegment(HitCleaningEHE, 'eheclean',
	inATWD='EHECalibratedATWD_Wave', inFADC = 'EHECalibratedFADC_Wave',
	If=which_split(split_name='InIceSplit') & (lambda f: ehe_wg(f))
	)

# tray.AddModule(TesterModule, 'tester', Streams=[icetray.I3Frame.DAQ])
# tray.Add("Dump")
# tray.Add("I3Writer", filename="134777_8912764_L2_standard.i3.zst")
# tray.Add("I3Writer", filename="134777_8912764_L2_FADC_True_ATWD_True.i3.zst")
# tray.Add("I3Writer", filename="134777_8912764_L2_FADC_True_ATWD_False.i3.zst")
tray.Add("I3Writer", filename="134777_8912764_L2_FADC_False_ATWD_True.i3.zst")
tray.Execute()