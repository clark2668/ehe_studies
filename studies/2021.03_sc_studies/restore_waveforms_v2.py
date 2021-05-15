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

from icecube import VHESelfVeto, hdfwriter
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
	help="exclude the atwd or not (True to exclude the atwd"
	)
args = parser.parse_args()
args.exclude_atwd = bool(distutils.util.strtobool(args.exclude_atwd))
print('Exclude ATWDs?: {}'.format(args.exclude_atwd))

# setup
#####################################
#####################################

tray = I3Tray()
scdata = '/misc/disk15/data/IceCube/RealData/86strings/standardcandle/2015/sc2/'

tray.AddModule("I3Reader", filenamelist=[scdata+'Level2_IC86.2014_data_Run00125920_0116_0_138_GCD.i3.gz', 
	args.input_file]) #GCD first

# we have to delete a bunch of stuff in order to get the Rehydration to run again successfully
tray.AddModule(drop_pframe, 'dropP')
tray.AddModule("Delete", 'deleter', 
	Keys=['InIceDSTPulses', 'IceTopDSTPulses', 'CleanInIceRawData', 'CleanIceTopRawData',
	'CalibratedWaveformRange', 'ReextractedInIcePulses', 'ReextractedInIcePulsesTimeRange',
	'ReextractedIceTopPulses', 'IceTopHLCPulseInfo', 'ReextractedIceTopPulses_SLC',
	'InIcePulses', 'IceTopPulses', 'RehydrateNInIcePFrames', 'NFramesIsDifferent', 'IceTopErrata',
	'CalibrationErrata', 'SaturationWindows', 'CalibratedSLC']
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
	excludeATWD=args.exclude_atwd
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



# we have to calculate a bunch of things
# we need HQtot and Portia, both "standard" and "magsix"
# for pulses unfolded using FADC only
# the "FADC"-only part is handled deep inside the L2 scripts

# HQtot first
#####################################
#####################################

# HQTot, standard mode
pulses='SplitInIcePulses'
tray.AddModule('HomogenizedQTot', 'qtot_total', 
	Pulses=pulses,
	Output='HomogenizedQTot'
	)

# HQtot, deep mag six
tray.AddModule(utils_pulses.CalcQTOt_DeepMagSix_module, 'hqtot_deepmagsix',
	pulses=pulses,
	name='HomogenizedQTot_DeepMagSix',
	Streams=[icetray.I3Frame.Physics],
	)

# Portia second
#####################################
#####################################

# bog-standard Portia

# Portia Best NPE, normal mode
tray.AddModule(utils_pulses.CalcPortiaCharge_module, 'portia_standard',
	excludeATWD=args.exclude_atwd, excludeFADC=False,
	name='PortiaEventSummarySRT',
	Streams=[icetray.I3Frame.Physics],
	)

# Portia Best NPE, deep mag six
tray.AddModule(utils_pulses.CalcPortiaCharge_DeepMagSix_module, 'portia_deepmagsix',
	excludeATWD=args.exclude_atwd, excludeFADC=False,
	name='PortiaEventSummarySRT_DeepMagSix',
	Streams=[icetray.I3Frame.Physics],
	)


# output
#####################################
#####################################

tray.AddSegment(hdfwriter.I3HDFWriter, 'hdf', 
	Output=f'{args.output_file}.hdf5', 
	Keys=['I3EventHeader', 'HomogenizedQTot', 'HomogenizedQTot_DeepMagSix',
	'PortiaEventSummarySRT', 'PortiaEventSummarySRT_DeepMagSix'], 
	SubEventStreams=['InIceSplit']
	)

# tray.Add("Dump")
tray.AddModule("I3Writer", "write",
	filename=f'{args.output_file}.i3.zst',
	Streams=[icetray.I3Frame.DAQ, icetray.I3Frame.Physics],
	)

tray.Execute(10)
