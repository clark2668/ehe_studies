#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/icetray-start
#METAPROJECT combo/V01-01-00

from icecube import icetray, dataio, dataclasses, common_variables, linefit
from icecube import phys_services
from I3Tray import I3Tray
from icecube import hdfwriter
import numpy as np


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


def calculate_closest_approach(frame, primary_name):
    primaryEvent = frame[primary_name]
    position = primaryEvent.pos
    direction = primaryEvent.dir
    origin = dataclasses.I3Position(0., 0., 0.)
    approach = (position + direction * (direction*(origin - position))).magnitude
    # just to validate with standard tool (only really works for numu)
    # approach_v2 = phys_services.I3Calculator.closest_approach_distance(primaryEvent,origin)
    # print("Closest Approach me {}, Calc {}".format(approach, approach_v2))
    frame["ClosestApproach"] = dataclasses.I3Double(approach)
        
tray.AddModule(calculate_closest_approach, 'calcClosest',
    primary_name = 'PrimaryEvent',
    Streams=[icetray.I3Frame.Physics]
    )
        

tray.AddSegment(hdfwriter.I3HDFWriter, 'hdf', 
    Output=f'{args.output_file}.hdf5', 
    Keys=[
        'I3EventHeader', 'CorsikaWeightMap', 'PolyplopiaPrimary', 'I3MCWeightDict',
        'CVMultiplicity', 'CVStatistics', 'Homogenized_QTot', 'EHELineFit',
        'LineFit', 'PrimaryEvent', 'I3JulietPrimaryParticle', 
        'PropagationMatrixNuE', 'PropagationMatrixNuMu', 'PropagationMatrixNuTau',
        'JulietWeightDict', 
        'EHEOpheliaParticleSRT_ImpLF', 'EHEOpheliaSRT_ImpLF', 'EHEPortiaEventSummarySRT',
        'ClosestApproach'
        
    ],
    SubEventStreams=['InIceSplit', 'Final']
    )

tray.Execute()