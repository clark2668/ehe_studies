#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/icetray-start
#METAPROJECT combo/V01-01-00

from icecube import icetray, dataio, dataclasses
from I3Tray import I3Tray
from icecube import VHESelfVeto, hdfwriter
import numpy as np

import sys
sys.path.append('/home/brianclark/IceCube/ehe_studies/utils/')

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-i", type=str,
    dest="input_file",required=True,
    help="full path to the input file",
    )
parser.add_argument("-o", type=str,
    dest="output_file",required=True,
    help='''full path to the output file, without file extention. 
            That is, provide "test" not "test.i3.bz2"''',
    )
parser.add_argument("-s", type=bool,
    dest="save_i3file",required=False, default=False,
    help='should save i3 file',
    )
parser.add_argument("-g", type=str,
    dest="gcd_file",required=False, default=None,
    help='gcd file',
    )
# /cvmfs/icecube.opensciencegrid.org/data/GCD/GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE_withScaledNoise.i3.gz
args = parser.parse_args()

filenamelist = []
if args.gcd_file is not None:
    filenamelist.append(args.gcd_file)
filenamelist.append(args.input_file)

tray = I3Tray()
tray.AddModule("I3Reader", filenamelist=filenamelist)

# cut on the high Q filter
from icecube.filterscripts import filter_globals
def highQfilter(frame):
    if frame.Stop == icetray.I3Frame.Physics and frame.Has('FilterMask'):
        if frame['FilterMask'].get(filter_globals.HighQFilter).condition_passed:
            return 1
        else:
            return 0
    else:
        return 0

tray.AddModule(highQfilter, 'highQ',
    Streams=[icetray.I3Frame.Physics])

do_fill_ratio = False
if do_fill_ratio:

    '''
    If we need to calculate the fill ratio ourselves
    '''

    suffix = '_L3'
    tray.AddModule('Delete', 'deleter', 
        Keys=['OfflinePulsesSLC', 'OfflinePulsesHLC', 'TWOfflinePulsesHLC']
        )

    # first, re-run cascade hit-cleaning
    from icecube.phys_services.which_split import which_split
    from icecube.filterscripts.offlineL2.level2_HitCleaning_Cascade import CascadeHitCleaning
    tray.AddSegment(CascadeHitCleaning, 'CascadeHitCleaning',
        If=which_split(split_name='InIceSplit')
        )

    # rerun some subset of the cascade reco
    import utils_reco
    from utils_reco import OfflineCascadeReco
    tray.AddSegment(OfflineCascadeReco, 'CascadeL2Reco',
        SRTPulses='SRTInIcePulses',
        Pulses='TWOfflinePulsesHLC',
        TopoPulses='OfflinePulsesHLC',
        If=which_split(split_name='InIceSplit'),
        suffix=suffix
        )

# def make_cut(frame):
#     hqtot = frame.Get('Homogenized_QTot').value
#     fillratio = frame.Get('CascadeFillRatio_L3').fill_ratio_from_mean
#     speed = frame.Get('LineFit').speed
#     if np.log10(hqtot) < 5:
#         return 0
#     if fillratio < 0.95:
#         return 0
#     if speed > 0.25:
#         return 0

#     return 1

# # make cuts
# tray.AddModule(make_cut, 'cut', Streams=[icetray.I3Frame.Physics])

def is_neutrino(pType):
       return(pType==dataclasses.I3Particle.ParticleType.NuE
              or pType==dataclasses.I3Particle.ParticleType.NuEBar
              or pType==dataclasses.I3Particle.ParticleType.NuMu
              or pType==dataclasses.I3Particle.ParticleType.NuMuBar
              or pType==dataclasses.I3Particle.ParticleType.NuTau
              or pType==dataclasses.I3Particle.ParticleType.NuTauBar)

def find_final_neutrino(frame, mctree_name):
        mcTree = frame[mctree_name]
        # print(mcTree)
        #try to figure out which neutrino is the primary
        primaryNeutrino=dataclasses.get_most_energetic_primary(mcTree)
        frame["PrimaryNeutrino"]=primaryNeutrino
        if(primaryNeutrino==None or not is_neutrino(primaryNeutrino.type)):
            return
        #walk down the tree to find the first daughter neutrino which is 'InIce'
        neutrino=primaryNeutrino
        while(neutrino.location_type!=dataclasses.I3Particle.LocationType.InIce):
            children=mcTree.get_daughters(neutrino)
            foundNext=False
            #take the first child which is a neutrino;
            #for in-Earth NC interactions it should be the only one anyway
            for child in children:
                if(is_neutrino(child.type)):
                    neutrino=child
                    foundNext=True
                    break
            if(not foundNext):
                print("did not find a daughter neutrino")
                return #bail out
        frame["InteractingNeutrino"]=neutrino
        stop_pos = neutrino.pos + neutrino.dir*neutrino.length
        frame["VertexPosition"] = stop_pos
        # if(neutrino==primaryNeutrino):
        #     print('Final selection is same as first, z pos {}'.format(stop_pos.z))
        # else:
        #     print("other, x/y/z {:e}, {:e}, {:e}".format(stop_pos.x, stop_pos.y, stop_pos.z))
        #     print(mcTree)


tray.AddModule(find_final_neutrino, 'vertex', mctree_name='I3MCTree_preMuonProp',
    Streams=[icetray.I3Frame.Physics]
    )

name = 'redo'
# redo_name = 'redo_LF'
input_pulses = 'SRTInIcePulses'

# redo linefit itself
from icecube import linefit
tray.AddSegment(linefit.simple, 'redo_LF', inputResponse=input_pulses,
    fitName='LineFit_'+name
    )

from utils_pulses import RedoLineFitPulseDebiasing
tray.AddSegment(RedoLineFitPulseDebiasing, name, input_pulses=input_pulses)

# get the fit quality of vanilla linefit
from utils_pulses import get_linefit_quality
tray.AddModule(get_linefit_quality, 'LFqual', 
    linefit_name='LineFit_'+name,
    linefit_params_name='LineFit_'+name+'Params',
    pulses_name=name+'_debiasedPulses',
    output_name='LineFit_'+name+'Quality'
    )

# then recalculate for an impact factor of 300
impact = 300
output_pulses1 = name+'_debiasedPulses' + '_CutFarAway_{}'.format(impact)
from utils_pulses import strip_pulses
tray.AddModule(strip_pulses, 'sp1', 
    track_name='LineFit_'+name, impact_parameter=impact,
    input_pulses_name=input_pulses, output_pulses_name=output_pulses1
    )

qualfit_output_name_1 = 'LineFit_'+name+'Quality_CutFarAway_{}'.format(impact)
tray.AddModule(get_linefit_quality, 'LFqual2', 
    linefit_name='LineFit_'+name,
    linefit_params_name='LineFit_'+name+'Params',
    pulses_name=output_pulses1,
    output_name=qualfit_output_name_1
    )

# and one more time for an impact factor of 120 (Ty's suggestion)
impact = 120
output_pulses2 = name+'_debiasedPulses' + '_CutFarAway_{}'.format(impact)
from utils_pulses import strip_pulses
tray.AddModule(strip_pulses, 'sp2', 
    track_name='LineFit_'+name, impact_parameter=impact,
    input_pulses_name=input_pulses, output_pulses_name=output_pulses2
    )

qualfit_output_name_2 = 'LineFit_'+name+'Quality_CutFarAway_{}'.format(impact)
tray.AddModule(get_linefit_quality, 'LFqual3', 
    linefit_name='LineFit_'+name,
    linefit_params_name='LineFit_'+name+'Params',
    pulses_name=output_pulses2,
    output_name=qualfit_output_name_2
    )

# we also need to count channels and strings hit, etc. (hit multiplicity)
from icecube.common_variables import hit_multiplicity
pulses = 'SplitInIcePulses'
hmoutput_name = 'HitMultiplicityValues'
tray.AddSegment(hit_multiplicity.I3HitMultiplicityCalculatorSegment, 'hm',
    PulseSeriesMapName=pulses,
    OutputI3HitMultiplicityValuesName = hmoutput_name,
    # BookIt = True
)

# and also need things like the COG (hit statistics)
from icecube.common_variables import hit_statistics
hsoutput_name = 'HitStatisticsValues'
tray.AddSegment(hit_statistics.I3HitStatisticsCalculatorSegment, 'hs',
    PulseSeriesMapName=pulses,
    OutputI3HitStatisticsValuesName=hsoutput_name
)



tray.AddSegment(hdfwriter.I3HDFWriter, 'hdf', 
    Output=f'{args.output_file}.hdf5', 
    Keys=['I3MCWeightDict', 'I3EventHeader', 'PolyplopiaPrimary', 'CorsikaWeightMap',
    'InteractingNeutrino', 'VertexPosition','CascadeFillRatio_L3', 'PrimaryNeutrino',
    'Homogenized_QTot', 'LineFit', 'LineFit_'+name,
    'EHEPortiaEventSummarySRT', 'EHEOpheliaParticleSRT_ImpLF', 'EHEOpheliaSRT_ImpLF',
    'LineFit_'+name+'Quality', qualfit_output_name_1, qualfit_output_name_2,
    hmoutput_name, hsoutput_name, 
    ], 
    SubEventStreams=['InIceSplit']
    )

if args.save_i3file:
    tray.AddModule("I3Writer", "write",
        filename=f'{args.output_file}.i3.zst',
        Streams=[icetray.I3Frame.TrayInfo, icetray.I3Frame.DAQ, icetray.I3Frame.Physics],
        DropOrphanStreams=[icetray.I3Frame.Calibration, icetray.I3Frame.DAQ]
        )

tray.Execute()
