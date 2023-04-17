#!/usr/bin/env python3

from icecube import icetray, dataio, dataclasses, common_variables, linefit, gulliver
from icecube import phys_services
from I3Tray import I3Tray
from icecube import hdfwriter, weighting_module

import sys
sys.path.append('/home/brian/IceCube/ehe/ehe_software/venv_ehe')

from eheanalysis import millipede

keeps = [
        'I3EventHeader', 'CorsikaWeightMap', 'PolyplopiaPrimary', 'I3MCWeightDict',
        'I3JulietPrimaryParticle', 'JulietWeightDict', 
        'PropagationMatrixNuE', 'PropagationMatrixNuMu', 'PropagationMatrixNuTau',
        'CVMultiplicity', 'CVStatistics', 'Homogenized_QTot', 
        'EHELineFit', 'LineFit', 'LineFitEHE',
        'EHEOpheliaParticleSRT_ImpLF', 'EHEOpheliaSRT_ImpLF', 'EHEPortiaEventSummarySRT',
        'EHE_SplineMPE', "EHE_Monopod"
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


args = parser.parse_args()

gcd_file = '/cvmfs/icecube.opensciencegrid.org/data/GCD/GeoCalibDetectorStatus_2020.Run134142.Pass2_V0.i3.gz'
convex_hull = millipede.get_convex_hull(gcd_file)
# polygon = millipede.get_extruded_polygon(gcd_file)

filenamelist = []
filenamelist.append(args.input_file)

tray = I3Tray()
tray.context['I3FileStager'] = dataio.get_stagers()
tray.AddModule("I3Reader", filenamelist=filenamelist)

def trash_tracks(frame, reco_name):
    if reco_name not in frame:
        return False
    else:
        return True
tray.Add(trash_tracks, 'trashTracks', reco_name='EHE_Monopod')

def find_primary(frame, mctree_name):
    try:
        mcTree = frame[mctree_name]
    except:
        mcTree = frame['I3MCTree_preMuonProp'] # get this if the other isn't available
    primaries = mcTree.primaries
    primaryEvent = primaries[0]
    frame["PrimaryEvent"] = primaryEvent

tray.AddModule(find_primary, 'findPrimary',
    mctree_name = 'I3MCTree',
    Streams=[icetray.I3Frame.Physics]
    )

# calculate EM equivalent visible energy (correct for hadronic)
tray.AddModule(millipede.calculate_em_equiv_dep_energy,
               mctree_name='I3MCTree',
               primary_name='PrimaryEvent',
               output_name="EMEquivVisDepE"
               )
keeps.extend(['EMEquivVisDepE'])

# straight up "visible" energy (no hadronic correction)
def calculate_depe(frame):
    dep_e = frame['NuEProperties']['visible_energy']
    frame['DepE'] = dataclasses.I3Double(dep_e)

tray.AddModule(calculate_depe, 'calcDepE')
keeps.extend(['DepE'])

# vertex distance from hull to check containment
# from ic3_labels.labels.utils import geometry
# def calculate_dist_to_hull(frame, polygon, reco_particle_name):
#     reco_particle = frame.Get(reco_particle_name)
#     # position = reco_particle.pos
#     # dist = geometry.distance_to_convex_hull(hull, position)
#     dist = polygon.GetDistanceToHull(reco_particle.pos, reco_particle.dir)
#     print(dist)

from ic3_labels.labels.utils import geometry
def determine_containment(frame, hull, reco_particle_name, result_name):
    part = frame.Get(reco_particle_name)
    is_inside = geometry.point_is_inside(hull,
                                            (part.pos.x, part.pos.y, part.pos.z))
    frame.Put(result_name, icetray.I3Bool(is_inside))

tray.AddModule(determine_containment, 'determineContainment', 
               hull=convex_hull,
               reco_particle_name='EHE_Monopod',
               result_name="EHE_Monopod_Containment"
               )
keeps.extend(['EHE_Monopod_Containment'])

tray.Add('Delete', 
    Keys=['PropagationMatrixNuE', 'PropagationMatrixNuMu', 'PropagationMatrixNuTau',
    'JulietWeightDict'
    ]
)

# okay, now juliet weighting for L2 files
tray.AddModule(weighting_module.I3PropagationMatrixFiller,
               output_name='PropagationMatrix',
               propagation_matrix_dir='/data/IceCube/JULIeT/propMtxZeus/BB/HDF5',
               split_flavors=True)

tray.AddSegment(hdfwriter.I3HDFWriter, 'hdf', 
    Output=f'{args.output_file}.hdf5', 
    Keys=keeps,
    SubEventStreams=['InIceSplit', 'Final']
    )

tray.Execute()
