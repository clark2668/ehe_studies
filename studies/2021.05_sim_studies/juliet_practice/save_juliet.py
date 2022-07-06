
import os
from os.path import expandvars

from icecube import icetray, dataio, dataclasses
from icecube import VHESelfVeto, hdfwriter
from I3Tray import I3Tray
from I3Tray import load

from icecube import juliet_interface
load("libweighting-module")
load("libc2j-icetray")

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
args = parser.parse_args()

filenamelist = []
if args.gcd_file is not None:
    filenamelist.append(args.gcd_file)
filenamelist.append(args.input_file)

# icetray.logging.set_level("TRACE")
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


# now to do juliet stuff
JAVA_CLASS_PATH = os.path.expandvars('$I3_SRC/juliet/java_lib/classes')
tray.AddService("I3JavaVMFactory", "java_vm",
                Options=[expandvars(f"-Djava.class.path={JAVA_CLASS_PATH}"),
                "-Xms512m", "-Xmx1024m"])

tray.AddModule('I3WeigherModuleJuliet', 'nu_weight',
                frameMCWeightName='JulietWeightDict',
                PrimaryParticleType=0, # 0: neutrino  1: atmosphericMuon
                NeutrinoModelList=np.arange(1, 38).tolist(),
                NDOMthres=0)

tray.AddModule('I3WeigherModuleJuliet', 'atmo_weight',
                frameMCWeightName='JulietWeightDictAtmo',
                frameCosmicRayEnergyWeightName='CosmicRayEnergyDist',
                PrimaryParticleType=0, # 0: neutrino  1: atmosphericMuon
                NDOMthres=0,
                Alpha=1.97, # The Elbert formula's parameter
                MuonEthreshold=1505, # The Elbert formula's parameter [GeV]
                CELFlag=False, # use the propagation matrix
                gzkCutOffFlag=False, # With NO GZK cutoff. But you will fill the cutoff case later anyway
                AtmFluxWeightName='ElbertModelIC22',
                widthOfLogEnergyBin=0.05 # width of the logE bin of the CR energy distribution
                )

tray.AddModule("I3PropagationMatrixFiller", "filler",
                framePropagationMatrixVectorName="PropagationMatrix",
                NDOMthres=0,
                splitFlavors=True,
                widthOfLogEnergyBin=0.05 # width of the logE bin
                                        # of the Propagation Matrix
)

outkeys = [ 'I3EventHeader', 'PolyplopiaPrimary', 'Homogenized_QTot', 'LineFit',
            'EHEPortiaEventSummarySRT', 'EHEOpheliaParticleSRT_ImpLF', 'EHEOpheliaSRT_ImpLF',
            'I3JulietPrimaryParticle', 'JulietWeightDict', 'JulietWeightDictAtmo',
            'PropagationMatrix', 'PropagationMatrixNuE', 'PropagationMatrixNuMu', 'PropagationMatrixNuTau'
        ]


# tray.AddSegment(hdfwriter.I3HDFWriter, 'hdf', 
#     Output=f'{args.output_file}.hdf5', 
#     Keys=outkeys, 
#     SubEventStreams=['InIceSplit']
#     )
from icecube import common_variables, gulliver, millipede, \
                    phys_services, recclasses, simclasses, \
                    linefit, spline_reco, cramer_rao, ddddr
from icecube.tableio import I3TableWriter
from icecube.hdfwriter import I3HDFTableService
out_file = f'{args.output_file}.hdf5'
service = I3HDFTableService(out_file)
tray.AddModule(I3TableWriter, 'writer_hd5',
    tableservice=[service],
    keys=outkeys,
    SubEventStreams=['InIceSplit']
)

if args.save_i3file:
    tray.AddModule("I3Writer", "write",
        filename=f'{args.output_file}.i3.zst',
        Streams=[icetray.I3Frame.TrayInfo, icetray.I3Frame.DAQ, icetray.I3Frame.Physics],
        DropOrphanStreams=[icetray.I3Frame.Calibration, icetray.I3Frame.DAQ]
        )

tray.Execute(10)
