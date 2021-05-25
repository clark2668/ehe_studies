#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/icetray-start
#METAPROJECT /home/brian/IceCube/modified_icetray/bld/

'''
This version of the waveform restorer is focused on FADC-only estimates
'''

from icecube import icetray
from icecube.icetray import I3Units
from I3Tray import I3Tray

from icecube.filterscripts.offlineL2.Rehydration import Rehydration
from icecube.filterscripts.offlineL2.level2_EHE_Calibration import EHECalibration
from icecube.filterscripts.offlineL2.level2_HitCleaning_EHE import HitCleaningEHE
from icecube.phys_services.which_split import which_split

from icecube import hdfwriter
from utils import utils_pulses
import argparse
import distutils.util

class drop_pframe(icetray.I3Module):
	def __init__(self, context):
		icetray.I3Module.__init__(self, context)
	# def Configure(self):
		# pass
	def Physics(self, frame):
		pass # drop P frames; leave everything else alone

def drop_nullsplit(frame):
	return frame['I3EventHeader'].sub_event_stream!='NullSplit'

parser = argparse.ArgumentParser()
parser.add_argument("-i", type=str,
	dest="input_file",
	help="full path to the input file",
	required=True
	)
parser.add_argument("-o", type=str, 
	dest="output_file",required=True,
	help="full path to the output file"
	)
parser.add_argument("-atwd", type=str, 
	dest="exclude_atwd",required=True,
	help="exclude the atwd or not (True to exclude the atwd)"
	)
parser.add_argument("-fadc", type=str, 
	dest="exclude_fadc",required=True,
	help="exclude the fadc or not (True to exclude the fadc)"
	)
args = parser.parse_args()
args.exclude_atwd = bool(distutils.util.strtobool(args.exclude_atwd))
args.exclude_fadc = bool(distutils.util.strtobool(args.exclude_fadc))
print('Exclude ATWDs?: {}'.format(args.exclude_atwd))
print('Exclude FADCs?: {}'.format(args.exclude_fadc))

# setup
#####################################
#####################################

tray = I3Tray()
scdata = '/misc/disk15/data/IceCube/RealData/86strings/standardcandle/2015/sc2/'

tray.AddModule("I3Reader", filenamelist=[scdata+'Level2_IC86.2014_data_Run00125920_0116_0_138_GCD.i3.gz', 
	args.input_file]) #GCD first

# ncomplete = 0
# def increment(frame):
# 	global ncomplete
# 	ncomplete +=1 
# 	print('Completed frame {}'.format(ncomplete))

# # quick cut on 2021-05-25 to expedite finding bright events to look at waveforms from
# def quick_cut(frame):
# 	keeper = False
# 	if frame.Has('HomogenizedQTot'):
# 		hqtot = frame.Get('HomogenizedQTot').value
# 		if hqtot > 1E3:
# 			keeper = True
# 	return keeper
# tray.AddModule(quick_cut, 'cut',
# 	Streams=[icetray.I3Frame.Physics]
# 	)
# tray.Add("I3OrphanQDropper") # nuke orphan Q frames now

# we have to delete a bunch of stuff in order to get the Rehydration to run again successfully
tray.AddModule(drop_pframe, 'dropP')
tray.AddModule("Delete", 'deleter', 
	Keys=['InIceDSTPulses', 'IceTopDSTPulses', 'CleanInIceRawData', 'CleanIceTopRawData',
	'CalibratedWaveformRange', 'ReextractedInIcePulses', 'ReextractedInIcePulsesTimeRange',
	'ReextractedIceTopPulses', 'IceTopHLCPulseInfo', 'ReextractedIceTopPulses_SLC',
	'InIcePulses', 'IceTopPulses', 'RehydrateNInIcePFrames', 'NFramesIsDifferent', 'IceTopErrata',
	'CalibrationErrata', 'SaturationWindows', 'CalibratedSLC', 'HLCOfflineCleanInIceRawDataWODC', 
	'EHECalibratedATWD_Wave', 'EHECalibratedFADC_Wave'
	]
	)

# rehydrate, recalibrate, and re-split
#####################################
#####################################

# do rehydration and recalibration
tray.AddSegment(Rehydration, 'rehydrator',
	dstfile=None,
	mc=False,
	doNotQify=True, # do not Qify, as we already have Q frames
	pass2=False,
	excludeATWD=args.exclude_atwd,
	excludeFADC=args.exclude_fadc
	)
tray.Add(drop_nullsplit)

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
    If = which_split(split_name='InIceSplit')
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
	If=which_split(split_name='InIceSplit') #& (lambda f: ehe_wg(f))
	)


keys = tray.Add(utils_pulses.CalcChargeStatistics, 
	excludeFADC=args.exclude_fadc,
	excludeATWD=args.exclude_atwd
	)

hdf_keys = keys + ['I3EventHeader']

# tray.Add(increment, Streams=[iceeray.I3Frame.DAQ])


# output
#####################################
#####################################

tray.AddSegment(hdfwriter.I3HDFWriter, 'hdf', 
	Output=f'{args.output_file}.hdf5', 
	Keys=hdf_keys, 
	SubEventStreams=['InIceSplit']
	)

# tray.Add("Dump")
tray.AddModule("I3Writer", "write",
	filename=f'{args.output_file}.i3.zst',
	Streams=[icetray.I3Frame.DAQ, icetray.I3Frame.Physics],
	)

tray.Execute()
