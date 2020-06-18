import h5py
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import ROOT
from ROOT import gStyle, gPad

gStyle.SetOptStat(0)

parser = argparse.ArgumentParser()
parser.add_argument("-f", type=str, nargs='+',
	dest="input_files",
	help="paths to input files (absolute or relative)")

args = parser.parse_args()
files = args.input_files

h1_zenith = ROOT.TH1D("h1_zenith","h1_zenith",20,-1,1)
h1_npe = ROOT.TH1D("h1_npe","h1_npe",20,4,6)

for file in files:
	print(file)
	file_in = h5py.File(file, "r")
	data = file_in['data']
	portia_npe = data['portia_npe']
	ophelia_zenith = data['ophelia_zenith']
	for event, this_portia_npe in enumerate(portia_npe):
		if(this_portia_npe>25000):
			h1_zenith.Fill(ophelia_zenith[event])
			h1_npe.Fill(np.log10(portia_npe[event]))


# h1 zenith
c_zenith = ROOT.TCanvas("c_zenith","c_zenith",1100,850)
c_zenith.cd()
h1_zenith.Draw("histe")
h1_zenith.GetXaxis().SetTitle("cos(#theta)")
h1_zenith.GetYaxis().SetTitle("Number of Events")
h1_zenith.SetTitle("")
gPad.SetLogy()
c_zenith.SaveAs('ehe_zenith_{:d}events.pdf'.format(int(h1_zenith.GetEntries())))

# h1 npe
c_npe = ROOT.TCanvas("c_npe","c_npe",1100,850)
c_npe.cd()
h1_npe.Draw("histe")
h1_npe.GetXaxis().SetTitle("log_{10}(NPE)")
h1_npe.GetYaxis().SetTitle("Number of Events")
h1_npe.SetTitle("")
gPad.SetLogy()
c_npe.SaveAs('ehe_npe_{:d}events.pdf'.format(int(h1_npe.GetEntries())))
