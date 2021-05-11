'''
Script to study why some DOMs have lots of portia charge by little hqtot
Probably deals with splits
'''

import matplotlib.pyplot as plt
import numpy as np
import copy

from icecube import icetray, dataio, dataclasses
from I3Tray import I3Tray
from utils import utils_pulses

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-i", type=str,
	dest="input_files",required=True,
	help="full path to the input file",
	)
args = parser.parse_args()

scdata = '/misc/disk15/data/IceCube/RealData/86strings/standardcandle/2015/sc2/'
gcdfile = 'Level2_IC86.2014_data_Run00125920_0116_0_138_GCD.i3.gz'

tray = I3Tray()
tray.AddModule("I3Reader", filenamelist=[scdata+gcdfile, args.input_files])

# tray.AddModule(utils_pulses.CalcPortiaCharge_module, 'portia')
# tray.AddModule(utils_pulses.CalcQTot_module, 'hqtot')
tray.AddModule(utils_pulses.Compare_Portia_QTot, 'compare')

tray.Execute()

