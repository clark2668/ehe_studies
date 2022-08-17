#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.1.0/icetray-start
#METAPROJECT: /home/mmeier/data/software/metaprojects/icetray_main/build

###!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.1.0/icetray-start
###METAPROJECT: combo/V00-00-03

###!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.1.0/icetray-start
###METAPROJECT: /home/mmeier/data/software/metaprojects/icetray_main/build

import time
import os
import copy
import numpy as np
from optparse import OptionParser

from icecube import icetray, dataio, dataclasses, common_variables, linefit, phys_services, hdfwriter
from I3Tray import I3Tray

from I3Tray import *
from icecube.filterscripts import filter_globals
from icecube import filter_tools

from icecube.icetray import I3Units

from icecube.filterscripts.offlineL2 import Globals
from icecube.filterscripts.offlineL2.Globals import (deepcore_wg, 
    muon_wg, wimp_wg, cascade_wg, 
    fss_wg, fss_wg_finiteReco, ehe_wg, ehe_wg_Qstream, monopole_wg)
from icecube.filterscripts.offlineL2.Rehydration import Rehydration, Dehydration
from icecube.filterscripts.offlineL2.level2_IceTop_CalibrateAndExtractPulses import CalibrateAndExtractIceTop
from icecube.filterscripts.offlineL2.level2_EHE_Calibration import EHECalibration
from icecube.filterscripts.offlineL2.level2_HitCleaning_IceTop import IceTopCoincTWCleaning
from icecube.filterscripts.offlineL2.level2_HitCleaning_DeepCore import DeepCoreHitCleaning
from icecube.filterscripts.offlineL2.level2_HitCleaning_WIMP import WimpHitCleaning
from icecube.filterscripts.offlineL2.level2_HitCleaning_Cascade import CascadeHitCleaning
from icecube.filterscripts.offlineL2.PhotonTables import InstallTables
from icecube.filterscripts.offlineL2.level2_Reconstruction_Muon import OfflineMuonReco
from icecube.filterscripts.offlineL2.level2_HitCleaning_EHE import HitCleaningEHE
from icecube.filterscripts.offlineL2.level2_Reconstruction_IceTop import ReconstructIceTop
from icecube.filterscripts.offlineL2.level2_Reconstruction_DeepCore import OfflineDeepCoreReco
from icecube.filterscripts.offlineL2.level2_Reconstruction_WIMP import WimpReco
from icecube.filterscripts.offlineL2.level2_Reconstruction_Cascade import OfflineCascadeReco
from icecube.filterscripts.offlineL2.level2_Reconstruction_SLOP import SLOPLevel2
from icecube.filterscripts.offlineL2.level2_Reconstruction_EHE import ReconstructionEHE
from icecube.filterscripts.offlineL2.level2_Reconstruction_Monopole import MonopoleL2
from icecube.filterscripts.offlineL2.ClipStartStop import ClipStartStop
from icecube.phys_services.which_split import which_split

parser = OptionParser()
parser.add_option("-i", "--input", dest="INPUT",
                  type=str,
                  help="Input i3 file")
gcd_default = os.path.join(
    # '/cvmfs/icecube.opensciencegrid.org/data/GCD',
    '/disk20/users/brian/IceCube/support_files',
    'GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE_withStdNoise.i3.gz')
parser.add_option("-g", "--gcd", dest="GCD",
                  type=str,
                  default=gcd_default,
                  help="GCD file for input i3 file")
parser.add_option("-o", "--outputfile", dest="OUTPUTFILE",
                  type=str,
                  help="Output i3 file")

def count_cal_errata(frame):
    if frame.Has('CalibrationErrata'):
        cal_errata = frame.Get('CalibrationErrata')
        len_cal_errata = len(cal_errata)
        frame['LenCalErrata'] = icetray.I3Int(len_cal_errata)
    else:
        frame['LenCalErrata'] = icetray.I3Int(0)
    header = frame.Get("I3EventHeader")
    runid = header.run_id
    evid = header.event_id
    thelen = frame.Get("LenCalErrata")
    print(f"Run {runid}, Evid {evid}, Cal Err: {thelen}")


def get_datasetid_and_runnumber(filename):
    basename = os.path.basename(filename) # this returns just the .i3 part
    basename = basename.replace('.i3.zst', '') # ditch the file ending
    basename = basename.replace('.i3.bz2', '')
    splitted = basename.split('.')
    dataset_id = splitted[2]
    run_number = splitted[3]
    return int(dataset_id), int(run_number)


def correct_run_ids(frame, input_file):
    header = copy.copy(frame['I3EventHeader'])
    del frame['I3EventHeader']
    # dataset_id = parse_datasetid_from_filename(input_file)
    # run_number = parse_runnumber_from_filename(input_file)
    dataset_id, run_number = get_datasetid_and_runnumber(input_file)
    new_run_id = int(dataset_id * 100000 + run_number)
    header.run_id = new_run_id
    frame['I3EventHeader'] = header
    return True

class drop_pframe(icetray.I3Module):
    def __init__(self, context):
        icetray.I3Module.__init__(self, context)
    # def Configure(self):
        # pass
    def Physics(self, frame):
        pass # drop P frames; leave everything else alone

def drop_nullsplit(frame):
    return frame['I3EventHeader'].sub_event_stream!='NullSplit'


def main(input_file, gcd_file, output_file):
    tray = I3Tray()
    # tray.context['I3FileStager'] = dataio.get_stagers()

    tray.AddModule('I3Reader', 'reader',
                   FilenameList=[gcd_file, input_file])
    print("Reading input file...", input_file)

    tray.AddModule(drop_nullsplit, 'dropNull')

    do_recalculate = False
    if do_recalculate:

        # clean up first
        
        tray.AddModule(drop_pframe, 'dropP')

        tray.AddModule("Delete", 'deleter', 
            Keys=['InIceDSTPulses', 'IceTopDSTPulses', 'CleanInIceRawData', 'CleanIceTopRawData',
            'CalibratedWaveformRange', 'ReextractedInIcePulses', 'ReextractedInIcePulsesTimeRange',
            'ReextractedIceTopPulses', 'IceTopHLCPulseInfo', 'ReextractedIceTopPulses_SLC',
            'InIcePulses', 'IceTopPulses', 'RehydrateNInIcePFrames', 'NFramesIsDifferent', 'IceTopErrata',
            'CalibrationErrata', 'SaturationWindows', 'CalibratedSLC'
            ]
            )
        
        # and then essentially do everything in the offline L2 script

        # do rehydration
        tray.AddSegment(Rehydration, 'rehydrator',
            dstfile=None,
            mc=False,
            doNotQify=True, # do not Qify, as we already have Q frames
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
            # If = which_split(split_name='InIceSplit') & (lambda f: (
            #                  deepcore_wg(f) or wimp_wg(f)    or
            #                  muon_wg(f)     or cascade_wg(f) or
            #                  ehe_wg(f)      or fss_wg(f) ))
            If = which_split(split_name='InIceSplit')
        )
        
        # EHE calibration
        tray.AddSegment(EHECalibration, 'ehecalib',
            inPulses='CleanInIceRawData',
            outATWD='EHECalibratedATWD_Wave',
            outFADC='EHECalibratedFADC_Wave',
            If=lambda f: ehe_wg_Qstream(f)
        )
        
        
        # cascade hit cleaning #
        tray.AddSegment(CascadeHitCleaning,'CascadeHitCleaning', 
            If=which_split(split_name='InIceSplit') & (lambda f: cascade_wg(f)),
        )
        
        # ehe hit cleaning #
        tray.AddSegment(HitCleaningEHE, 'eheclean',
            inATWD='EHECalibratedATWD_Wave', inFADC = 'EHECalibratedFADC_Wave',
            If=which_split(split_name='InIceSplit') & (lambda f: ehe_wg(f))
        )
        
        # load tables #
        # photonicsdir = '/disk20/users/brian/IceCube/support_files/photon-tables'
        photonicsdir = None
        if photonicsdir is not None:
            tray.AddSegment(InstallTables, 'InstallPhotonTables',
                PhotonicsDir=photonicsdir
            )
        else:
            tray.AddSegment(InstallTables, 'InstallPhotonTables')

        # muon, cascade, wimp, fss #
        tray.AddSegment(OfflineMuonReco, 'OfflineMuonRecoSLC',
            Pulses = "SRTInIcePulses",
            If = which_split(split_name='InIceSplit') & (lambda f: (
                        muon_wg(f) or cascade_wg(f) or
                        wimp_wg(f) or fss_wg(f))),
            suffix = "", #null? copied from level2_globals supplied
        )

        # cascade #
        tray.AddSegment(OfflineCascadeReco,'CascadeL2Reco',
            SRTPulses='SRTInIcePulses',
            Pulses='TWOfflinePulsesHLC',
            TopoPulses = 'OfflinePulsesHLC',
            If=which_split(split_name='InIceSplit') & (lambda f: cascade_wg(f)),
            suffix='_L2'
        )
        
        # ehe #
        tray.AddSegment(ReconstructionEHE, 'ehereco',
            Pulses='EHETWCInIcePulsesSRT',
            suffix='EHE', LineFit = 'LineFit',
            SPEFitSingle='SPEFitSingle', SPEFit = 'SPEFit12',
            N_iter=12,
            If=which_split(split_name='InIceSplit') & (lambda f: ehe_wg(f))
        )

    tray.AddModule(drop_nullsplit, 'dropNull2')

    tray.AddModule(correct_run_ids, 'correct_run_ids',
                input_file=input_file,
                Streams=[icetray.I3Frame.DAQ]
                )

    tray.AddModule(count_cal_errata, 'count_cal_errata',
                   Streams=[icetray.I3Frame.Physics]
                   )
    
    if do_recalculate:
        # when we're just reading out the old one, no need to rewrite the i3 files
        tray.AddModule('I3Writer', 'writer',
            # DropOrphanStreams=[icetray.I3Frame.DAQ],
            Streams=[icetray.I3Frame.TrayInfo, icetray.I3Frame.DAQ, icetray.I3Frame.Physics],
            Filename=output_file)

    hdf_name = output_file.replace('.i3.bz2', '.hdf5')
    hdf_name = output_file.replace('.i3.zst', '.hdf5')
    tray.AddSegment(hdfwriter.I3HDFWriter, 'hdf',
                    Output=hdf_name,
                    Keys=['I3EventHeader', 'LenCalErrata',
                          'PolyplopiaPrimary', 'I3MCWeightDict', 'Homogenized_QTot'
                          ],
                    SubEventStreams=['InIceSplit', 'Final']
                    )
    

    tray.AddModule('TrashCan', 'thecan')

    tray.Execute()
    tray.Finish()

    del tray


if __name__ == '__main__':
    (options, args) = parser.parse_args()
    input_file = options.INPUT
    output_file = options.OUTPUTFILE
    gcd_file = options.GCD

    start_time = time.asctime()
    print('Started:', start_time)

    main(input_file, gcd_file, output_file)

    stop_time = time.asctime()

    print('Started:', start_time)
    print('Ended:', stop_time)

