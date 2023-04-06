#!/usr/bin/env python3

from multiprocessing.sharedctypes import Value
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

filenamelist = []
filenamelist.append(args.input_file)

tray = I3Tray()
tray.context['I3FileStager'] = dataio.get_stagers()
tray.AddModule("I3Reader", filenamelist=filenamelist)

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

tray.AddModule(millipede.calculate_em_equiv_dep_energy,
               mctree_name='I3MCTree',
               primary_name='PrimaryEvent',
               output_name="EMEquivVisDepE"
               )
keeps.extend(['EMEquivVisDepE'])

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
