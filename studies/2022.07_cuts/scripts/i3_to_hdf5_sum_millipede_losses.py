#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/icetray-start
#METAPROJECT /home/brian/IceCube/ehe/software/build_icetray

from multiprocessing.sharedctypes import Value
from icecube import icetray, dataio, dataclasses, common_variables, linefit
from icecube import phys_services
from I3Tray import I3Tray
from icecube import hdfwriter
import numpy as np

import sys
sys.path.append('/home/brian/IceCube/ehe/max_tools/build')
from eheanalysis import millipede

gcd_file = '/cvmfs/icecube.opensciencegrid.org/data/GCD/GeoCalibDetectorStatus_2020.Run134142.Pass2_V0.i3.gz'
convex_hull = millipede.get_convex_hull(gcd_file)


keeps = [
        'I3EventHeader', 'CorsikaWeightMap', 'PolyplopiaPrimary', 'I3MCWeightDict',
        'CVMultiplicity', 'CVStatistics', 'Homogenized_QTot', 'EHELineFit',
        'LineFit', 'PrimaryEvent', 'I3JulietPrimaryParticle', 
        'PropagationMatrixNuE', 'PropagationMatrixNuMu', 'PropagationMatrixNuTau',
        'JulietWeightDict', 
        'EHEOpheliaParticleSRT_ImpLF', 'EHEOpheliaSRT_ImpLF', 'EHEPortiaEventSummarySRT',
        'ClosestApproach', 'LenCalErrata', 'LenSatWindows', 'EHE_SplineMPE'
    ]

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
parser.add_argument("-c", type=str,
    dest="corsika_selection", required=False,
    help="If you are processing Max's high energy corsika, do you want to select for proton or iron?"
    )


args = parser.parse_args()

select_this_cor_species = None
if args.corsika_selection:
    cor_sel = args.corsika_selection.lower() # reduce to lower case
    if not cor_sel in ["p", "fe"]:
        raise ValueError(f"Corsika species {cor_sel} not supported")
    if cor_sel == "p":
        select_this_cor_species = dataclasses.I3Particle.PPlus
    elif cor_sel == "fe":
        select_this_cor_species = dataclasses.I3Particle.Fe56Nucleus
    print(f"Selecting for corsika species {select_this_cor_species}")
    

filenamelist = []
filenamelist.append(args.input_file)

tray = I3Tray()
tray.context['I3FileStager'] = dataio.get_stagers()
tray.AddModule("I3Reader", filenamelist=filenamelist)

def filter_corsika_primary(frame, select_this_type):
    
    proton = dataclasses.I3Particle.PPlus
    iron = dataclasses.I3Particle.Fe56Nucleus
        
    if frame.Has("CorsikaWeightMap"):
        weight_map = frame.Get("CorsikaWeightMap")
        primary_type = int(weight_map['PrimaryType'])
        if(primary_type == select_this_type ):
            print("Keep {}".format(primary_type))
            return 1
        else:
            print("Trash {}".format(primary_type))
            return 0
    else:
        return 0

if select_this_cor_species is not None:
    tray.AddModule(filter_corsika_primary, 'filt_cor',
        select_this_type = select_this_cor_species,
        Streams=[icetray.I3Frame.DAQ, icetray.I3Frame.Physics]
        )

# high Q filter (must have >1000 PE)
def highQfilter(frame):
    try:
        charge = frame.Get("Homogenized_QTot")
        if charge >= 1E3:
            return 1
        else:
            return 0
    except:
        return 0

tray.AddModule(highQfilter, 'highQ',
    Streams=[icetray.I3Frame.Physics],
    )

def find_primary(frame, mctree_name):
    try:
        mcTree = frame[mctree_name]
    except:
        mcTree = frame['I3MCTree_preMuonProp'] # get this if the other isn't available
    primaries = mcTree.primaries
    primaryEvent = primaries[0]
    frame["PrimaryEvent"] = primaryEvent

tray.AddModule(find_primary, 'findNeutrino',
    mctree_name = 'I3MCTree',
    Streams=[icetray.I3Frame.Physics]
    )
millipede_name='EHEMuMillipede_SplineMPEseed'
tray.AddModule(millipede.SumMillipedeEnergyUnfolding,
    millipede_name=f'{millipede_name}',
    output_name=f'{millipede_name}_sum'
    )
keeps.append(f'{millipede_name}_sum')
tray.AddModule(millipede.SumMillipedeEnergyUnfolding,
    millipede_name=f'{millipede_name}',
    output_name=f'{millipede_name}_sum_contained',
    contained_only=True,
    convex_hull = convex_hull
    )
keeps.append(f'{millipede_name}_sum_contained')

tray.AddModule('I3Writer', 'writer',
    # DropOrphanStreams=[icetray.I3Frame.DAQ],
    Streams=[icetray.I3Frame.TrayInfo, icetray.I3Frame.DAQ, 
             icetray.I3Frame.Physics, icetray.I3Frame.Simulation],
    Filename=f'{args.output_file}.i3.zst')


tray.AddSegment(hdfwriter.I3HDFWriter, 'hdf', 
    Output=f'{args.output_file}.hdf5', 
    Keys=keeps,
    SubEventStreams=['InIceSplit', 'Final']
    )

tray.Execute()