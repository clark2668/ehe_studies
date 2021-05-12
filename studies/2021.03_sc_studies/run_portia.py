# this redoes *just* the portia bit of calibration / calculation
# if you need everything, or don't have the intermediate data products,
# better to use restore_waveforms.py

from icecube import icetray
from icecube.icetray import I3Units
from I3Tray import I3Tray

from icecube.filterscripts.offlineL2.Rehydration import Rehydration
from icecube.filterscripts.offlineL2.level2_EHE_Calibration import EHECalibration
from icecube.filterscripts.offlineL2.level2_HitCleaning_EHE import HitCleaningEHE
from icecube.filterscripts.offlineL2.Globals import ehe_wg, ehe_wg_Qstream

from icecube.phys_services.which_split import which_split
# icetray.set_log_level(icetray.I3LogLevel.LOG_INFO)

from icecube import VHESelfVeto, hdfwriter
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", type=str,
	dest="input_files",
	help="full path to the input file",
	required=True
	)
parser.add_argument("-o", type=str, 
	dest="output_dir",required=True,
	help="directory where the output should be written to"
	)
parser.add_argument("-y", type=int,
	dest="year", default=2015,
	help="Which year of standard candle, e.g. 2015"
	)
parser.add_argument("-c", type=int,
	dest="candle", default=2,
	help="Which candle, e.g. 2"
	)
parser.add_argument("-f", type=int,
	dest="filter_setting", required=True,
	help="Which standard candle filter setting, e.g. 0, 1, 2, ..."
	)
args = parser.parse_args()

tray = I3Tray()
scdata = '/misc/disk15/data/IceCube/RealData/86strings/standardcandle/2015/sc2/'
# tray.AddModule("I3Reader", filenamelist=[scdata+'Level2_IC86.2014_data_Run00125920_0116_0_138_GCD.i3.gz', 
	# scdata+'Level2_Run00125920_Subrun00000000_00000000.i3.bz2']) #GCD first

tray.AddModule("I3Reader", filenamelist=[scdata+'Level2_IC86.2014_data_Run00125920_0116_0_138_GCD.i3.gz', 
	args.input_files]) #GCD first

# we have to delete a bunch of stuff in order to get the Rehydration to run again successfully
tray.AddModule("Delete", 'deleter', 
	Keys=['HLCOfflineCleanInIceRawDataWODC', 'CleanIceTopRawData_EHE', 
	'EHEHLCCalibratedWaveforms', 'EHECalibratedATWD_Wave', 'EHECalibratedFADC_Wave',
	'CalibratedSLC', 'HLCOfflineCleanInIceRawData', 'SLCOfflineCleanInIceRawData',
	'splittedDOMMap', 'EHEPortiaEventSummary', 'LargestOMKey', 'splittedDOMMapSRT', 
	'EHEPortiaEventSummarySRT', 'EHEFADCPortiaPulseSRT', 'EHEFADCPulseSeriesSRT',
	'EHEATWDPortiaPulseSRT', 'EHEATWDPulseSeriesSRT', 'EHEBestPortiaPulseSRT',
	'SRTInIcePulses_WODC', 'SRTInIcePulses_WODCCleanedKeys', 'EHEInIcePulsesSRT',
	'EHETWCInIcePulsesSRT', 'HuberFit', 'EHEdebiased_BestPortiaPulseSRT_CleanDelay']
	)

# EHE Calibration
tray.AddSegment(EHECalibration, 'ehecalib',
	inPulses='CleanInIceRawData',
	outATWD='EHECalibratedATWD_Wave',
	outFADC='EHECalibratedFADC_Wave',
	# If=lambda f: ehe_wg_Qstream(f)
	)

# EHE Hit Cleaning
tray.AddSegment(HitCleaningEHE, 'eheclean',
	inATWD='EHECalibratedATWD_Wave', inFADC = 'EHECalibratedFADC_Wave',
	# inATWD='CalibratedWaveforms_ATWD', inFADC = 'CalibratedWaveforms_FADC',
	If=which_split(split_name='InIceSplit') #& (lambda f: ehe_wg(f))
	)
tray.AddModule("Delete", 'deleter2', 
	Keys=['HomogenizedQTot']
	)

pulses='SplitInIcePulses'
tray.AddModule('HomogenizedQTot', 'qtot_total', 
	Pulses=pulses,
	Output='HomogenizedQTot'
	)

# tray.Add("Dump")
tray.Add("I3Writer", 
	filename="{}/y{}_c{}_f{}_waves.i3.bz2".format(args.output_dir, args.year, args.candle, args.filter_setting),
	Streams=[icetray.I3Frame.DAQ, icetray.I3Frame.Physics]
	)

tray.Execute()
