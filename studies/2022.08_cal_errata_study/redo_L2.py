#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.1.0/icetray-start
#METAPROJECT: /home/mmeier/data/software/metaprojects/icetray_main/build

import time
import os
import copy
from optparse import OptionParser

from I3Tray import *
from icecube import dataio, phys_services
from icecube import icetray, dataclasses
from icecube.filterscripts.all_filters import OnlineFilter
from icecube.filterscripts import filter_globals
from icecube import filter_tools
from icecube.phys_services.which_split import which_split
from icecube.filterscripts.offlineL2.level2_all_filters import OfflineFilter

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


# def parse_runnumber_from_filename(filename):
#     basename = os.path.basename(filename)
#     basename = basename.replace('.i3.zst', '')
#     basename = basename.replace('.i3.bz2', '')
#     run_id = int(basename.split('_')[-1])
#     return run_id


# def parse_datasetid_from_filename(filename):
#     dirname = os.path.dirname(filename)
#     splitted = dirname.split('/')
#     flavor = splitted[-4]
#     energy_regime = splitted[-3]
#     flavor_to_int = {
#         'mu': 13,
#         'tau': 15,
#         'nue': 12,
#         'numu': 14,
#         'nutau': 15
#     }
#     energy_regime_to_int = {
#         'high_energy': 0,
#         'very_high_energy': 1
#     }
#     dataset_id = (flavor_to_int[flavor] * 1000 +
#         energy_regime_to_int[energy_regime])
#     return dataset_id

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
    print(frame['I3EventHeader'].sub_event_stream)
    print(frame['I3EventHeader'].sub_event_stream!='NullSplit')
    print("-----")
    # if frame['I3EventHeader'].sub_event_stream!='NullSplit'
    return frame['I3EventHeader'].sub_event_stream!='NullSplit'


def main(input_file, gcd_file, output_file):
    tray = I3Tray()

    tray.AddModule('I3Reader', 'reader',
                   FilenameList=[gcd_file, input_file])
    print("Reading input file...", input_file)

    # if 'I3_DATA' in os.environ:
    #     table_dir = os.path.expandvars("$I3_DATA/photon-tables/")
    # elif os.path.isdir('/cvmfs/icecube.opensciencegrid.org/data'):
    #     table_dir = "/cvmfs/icecube.opensciencegrid.org/data/photon-tables/"
    if os.path.isdir('/disk20/users/brian/IceCube/support_files'):
        table_dir = "/disk20/users/brian/IceCube/support_files/photon-tables/"
    else:
        raise Exception(
            'Cannot find I3_DATA or cvmfs, no photon tables available')

    spline_dir = os.path.join(table_dir, "splines")
    spline_reco_amplitude_table = os.path.join(
        spline_dir,
        'InfBareMu_mie_abs_z20a10_V2.fits')
    spline_reco_timing_table = os.path.join(
        spline_dir,
        'InfBareMu_mie_prob_z20a10_V2.fits')

    # clean up first
    
    tray.AddModule(drop_pframe, 'dropP')

    tray.AddModule(drop_nullsplit, 'dropNull')

    tray.AddModule("Delete", 'deleter', 
        Keys=['InIceDSTPulses', 'IceTopDSTPulses', 'CleanInIceRawData', 'CleanIceTopRawData',
        'CalibratedWaveformRange', 'ReextractedInIcePulses', 'ReextractedInIcePulsesTimeRange',
        'ReextractedIceTopPulses', 'IceTopHLCPulseInfo', 'ReextractedIceTopPulses_SLC',
        'InIcePulses', 'IceTopPulses', 'RehydrateNInIcePFrames', 'NFramesIsDifferent', 'IceTopErrata',
        'CalibrationErrata', 'SaturationWindows', 'CalibratedSLC', 'UncleanedInIcePulsesTimeRange',
        'I3SuperDST'
        ]
        )
    
    tray.AddSegment(OnlineFilter, "OnlineFilter",
                    sdstarchive=False, # our input data is SDST archive data
                    SplineRecoAmplitudeTable=spline_reco_amplitude_table,
                    SplineRecoTimingTable=spline_reco_timing_table,
                    simulation=True,
                    decode=False,
                    slop_split_enabled=False,
                    vemcal_enabled=False,
                    gfu_enabled=False,
                    needs_wavedeform_spe_corr=False,
                    alert_followup=False)

    # # make random service
    # seed = os.getpid()
    # filter_mask_randoms = phys_services.I3GSLRandomService(seed)

    # # override MinBias Prescale
    # filterconfigs = filter_globals.filter_pairs + filter_globals.sdst_pairs
    # print(filterconfigs)

    # # Generate filter Masks for all P frames
    # tray.AddModule(filter_tools.FilterMaskMaker, "MakeFilterMasks",
    #                OutputMaskName=filter_globals.filter_mask,
    #                FilterConfigs=filterconfigs,
    #                RandomService=filter_mask_randoms)

    # # Merge the FilterMasks
    # tray.AddModule("OrPframeFilterMasks", "make_q_filtermask",
    #                InputName=filter_globals.filter_mask,
    #                OutputName=filter_globals.qfilter_mask)

    # #Q+P frame specific keep module needs to go first, as KeepFromSubstram
    # #will rename things, let's rename post keep.
    # def is_Q(frame):
    #     return frame.Stop==frame.DAQ

    # gcd_keeps = ['I3Geometry',
    #              'I3Calibration',
    #              'I3DetectorStatus',
    #              'BadDomsList',
    #              'BadDomsListSLC',
    #              'BadDomsListHLC']

    # simulation_keeps = [
    #     'BackgroundI3MCTree',
    #     'BackgroundI3MCTreePEcounts',
    #     'BackgroundI3MCPESeriesMap',
    #     'BackgroundI3MCTree_preMuonProp',
    #     'BackgroundMMCTrackList',
    #     'BeaconLaunches',
    #     'CorsikaInteractionHeight',
    #     'CorsikaWeightMap',
    #     'EventProperties',
    #     'GenerationSpec',
    #     'I3LinearizedMCTree',
    #     'I3MCTree',
    #     'I3MCTreePEcounts',
    #     'I3MCTree_preMuonProp',
    #     'I3MCPESeriesMap',
    #     'I3MCPulseSeriesMap',
    #     'I3MCPulseSeriesMapParticleIDMap',
    #     'I3MCWeightDict',
    #     'LeptonInjectorProperties',
    #     'MCHitSeriesMap',
    #     'MCPrimary',
    #     'MCTrack',
    #     'MCPrimaryInfo',
    #     'MMCTrackList',
    # 'PolyplopiaCount',
    #     'PolyplopiaInfo',
    #     'PolyplopiaPrimary',
    #     'RNGState',
    #     'SignalI3MCPEs',
    #     'SimTrimmer', # for SimTrimmer flag
    #     'TimeShift', # the time shift amount
    #     'WIMP_params', # Wimp-sim
    #     'EventsPerFile',
    #     'I3JulietParamsTree',
    #     'I3JulietPrimaryParams',
    #     'I3JulietPrimaryParticle',
    #     'JulietRNGState'
    # ]

    # keep_before_merge = filter_globals.q_frame_keeps + [
    #     'InIceDSTPulses', # keep DST pulse masks
    #     'IceTopDSTPulses',
    #     'CalibratedWaveformRange', # keep calibration info
    #     'UncleanedInIcePulsesTimeRange',
    #     'SplitUncleanedInIcePulses',
    #     'SplitUncleanedInIcePulsesTimeRange',
    #     'SplitUncleanedInIceDSTPulsesTimeRange',
    #     'CalibrationErrata',
    #     'SaturationWindows',
    #     'InIceRawData', # keep raw data for now
    #     'IceTopRawData'] + simulation_keeps + gcd_keeps

    # tray.AddModule("Keep", "keep_before_merge",
    #                keys=keep_before_merge,
    #                If=is_Q)

    # # second set of prekeeps,
    # # conditional on filter content, based on newly created Qfiltermask
    # # Determine if we should apply harsh keep for events
    # # that failed to pass any filter
    # # Note: excluding the sdst_streams entries
    # tray.AddModule("I3IcePickModule<FilterMaskFilter>","filterMaskCheckAll",
    #                FilterNameList=filter_globals.filter_streams,
    #                FilterResultName=filter_globals.qfilter_mask,
    #                DecisionName="PassedAnyFilter",
    #                DiscardEvents=False,
    #                Streams=[icetray.I3Frame.DAQ])

    # def do_save_just_superdst(frame):
    #     if frame.Has("PassedAnyFilter"):
    #         if not frame["PassedAnyFilter"].value:
    #             # Event failed to pass any filter.
    #             return True
    #         else:
    #             # Event passed some filter
    #             return False
    #     else:
    #         print("Failed to find key frame Bool!!")
    #         return False

    # keep_only_superdsts = filter_globals.keep_nofilterpass + [
    #     'PassedAnyFilter',
    #     'InIceDSTPulses',
    #     'IceTopDSTPulses',
    #     'SplitUncleanedInIcePulses',
    #     'SplitUncleanedInIcePulsesTimeRange',
    #     'SplitUncleanedInIceDSTPulsesTimeRange',
    #     'RNGState'] + simulation_keeps + gcd_keeps

    # tray.AddModule("Keep", "KeepOnlySuperDSTs",
    #                keys=keep_only_superdsts,
    #                If=do_save_just_superdst)

    # # Now clean up the events that not even the SuperDST filters passed on.
    # tray.AddModule("I3IcePickModule<FilterMaskFilter>","filterMaskCheckSDST",
    #                FilterNameList=filter_globals.sdst_streams,
    #                FilterResultName=filter_globals.qfilter_mask,
    #                DecisionName="PassedKeepSuperDSTOnly",
    #                DiscardEvents=False,
    #                Streams=[icetray.I3Frame.DAQ])

    # def dont_save_superdst(frame):
    #     if frame.Has("PassedKeepSuperDSTOnly") and frame.Has("PassedAnyFilter"):
    #         if frame["PassedAnyFilter"].value:
    #             # These passed a regular filter, keep em
    #             return False
    #         elif not frame["PassedKeepSuperDSTOnly"].value:
    #             # Event failed to pass SDST filter.
    #             return True
    #         else:
    #             # Event passed some SDST filter
    #             return False
    #     else:
    #         print("Failed to find key frame Bool!!")
    #         return False

    # tray.AddModule("Keep", "KeepOnlyDSTs",
    #                keys=filter_globals.keep_dst_only +
    #                     ["PassedAnyFilter",
    #                      "PassedKeepSuperDSTOnly",
    #                      filter_globals.eventheader] +
    #                     gcd_keeps,
    #                If=dont_save_superdst)

    # # Frames should now contain only what is needed.
    # # now flatten, write/send to server
    # # Squish P frames back to single Q frame, one for each split:
    # tray.AddModule("KeepFromSubstream", "null_stream",
    #                StreamName=filter_globals.NullSplitter,
    #                KeepKeys=filter_globals.null_split_keeps)

    # in_ice_keeps = filter_globals.inice_split_keeps + filter_globals.onlinel2filter_keeps
    # in_ice_keeps = in_ice_keeps + [
    #     'I3EventHeader',
    #     'SplitUncleanedInIcePulses',
    #     'SplitUncleanedInIcePulsesTimeRange',
    #     'TriggerSplitterLaunchWindow',
    #     'I3TriggerHierarchy',
    #     'GCFilter_GCFilterMJD'] + gcd_keeps

    # tray.AddModule("Keep", "inice_keeps",
    #                keys=in_ice_keeps,
    #                If=which_split(split_name=filter_globals.InIceSplitter))

    # tray.AddModule("KeepFromSubstream", "icetop_split_stream",
    #                StreamName=filter_globals.IceTopSplitter,
    #                KeepKeys=filter_globals.icetop_split_keeps)

    # # Apply small keep list (SuperDST/SmallTrig/DST/FilterMask
    # # for non-filter passers
    # # Remove I3DAQData object for events not
    # # passing one of the 'filters_keeping_allraw'
    # tray.AddModule("I3IcePickModule<FilterMaskFilter>", "filterMaskCheck",
    #                FilterNameList=filter_globals.filters_keeping_allraw,
    #                FilterResultName=filter_globals.qfilter_mask,
    #                DecisionName="PassedConventional",
    #                DiscardEvents=False,
    #                Streams=[icetray.I3Frame.DAQ])

    # ## Clean out the Raw Data when not passing conventional filter
    # def I3RawDataCleaner(frame):
    #     if not (('PassedConventional' in frame and
    #              frame['PassedConventional'].value == True) or
    #             ('SimTrimmer' in frame and
    #              frame['SimTrimmer'].value == True)
    #            ):
    #         frame.Delete('InIceRawData')
    #         frame.Delete('IceTopRawData')

    # tray.AddModule(I3RawDataCleaner, "CleanErrataForConventional",
    #                Streams=[icetray.I3Frame.DAQ])

    # # DO L2 PROCESSING HERE
    # tray.AddSegment(OfflineFilter, "OfflineFilter",
    #     dstfile=False,
    #     mc=True,
    #     doNotQify=True,
    #     photonicsdir=table_dir)

    # def streamcut(frame):
    #     return frame["I3EventHeader"].sub_event_stream!='IceTopSplit' and \
    #         frame["I3EventHeader"].sub_event_stream!='NullSplit' and \
    #         frame["I3EventHeader"].sub_event_stream!='null_split'

    # tray.Add(streamcut)

    # tray.Add('Delete', 'final_delete',
    #          Keys=['ClusterCleaningExcludedTanks',
    #                'IceTopDSTPulses',
    #                'IceTopPulses',
    #                'CleanIceTopRawData',
    #                'IceTopRawData',
    #                'OfflineIceTopHLCTankPulses',
    #                'OfflineIceTopHLCVEMPulses',
    #                'OfflineIceTopSLCVEMPulses',
    #                'SimTrimmer',
    #                'TankPulseMergerExcludedTanks',
    #                'I3MCPESeriesMap',
    #                'I3MCPESeriesMap_081',
    #                'I3MCPESeriesMap_090',
    #                'I3MCPESeriesMap_095',
    #                'I3MCPESeriesMap_099',
    #                'I3MCPESeriesMap_108',
    #                'I3MCPESeriesMap_117'])
    #                #'InIceRawData'

    tray.AddModule(correct_run_ids, 'correct_run_ids',
                   input_file=input_file)

    tray.AddModule('I3Writer', 'writer',
        DropOrphanStreams=[icetray.I3Frame.DAQ],
        Streams=[icetray.I3Frame.DAQ, icetray.I3Frame.Physics],
        Filename=output_file)

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

