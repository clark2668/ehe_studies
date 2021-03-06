import h5py
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import ROOT
from ROOT import gStyle, gPad
from ROOT import kRed

gStyle.SetOptStat(0)

parser = argparse.ArgumentParser()
parser.add_argument("-f", type=str, nargs='+',
	dest="input_files",
	help="paths to input files (absolute or relative)")

args = parser.parse_args()
files = args.input_files

h2_npe_vs_hqt = ROOT.TH2D("h2","h2",400,2,6,400,2,6)
h2_npe_vs_hqt_highq = ROOT.TH2D("h2_npe_vs_hqt_highq","h2_npe_vs_hqt_highq",250,3.5,6,250,3.5,6)
h1_diff_npe_hqt = ROOT.TH1D("h1_diff_npe_hqt","h1_diff_npe_hqt",40,-2,2)

for file in files:
	print(file)
	file_in = h5py.File(file, "r")
	data = file_in['data']
	portia_npe = data['portia_npe']
	homogenized_qtot = data['homogenized_qtot']
	for event, this_portia_npe in enumerate(portia_npe):
		h2_npe_vs_hqt.Fill(np.log10(this_portia_npe), np.log10(homogenized_qtot[event]))
		if(np.log10(this_portia_npe)>3.5):
			h2_npe_vs_hqt_highq.Fill(np.log10(this_portia_npe), np.log10(homogenized_qtot[event]))
			if(np.log10(this_portia_npe)>4):
				toplot =  (this_portia_npe-homogenized_qtot[event])/this_portia_npe
				h1_diff_npe_hqt.Fill(toplot)

# draw a 1-1 line
line = ROOT.TLine(2,2,6,6)

# first th2d plot
c_npe_vs_hqt = ROOT.TCanvas("c_npe_vs_hqt","c_npe_vs_hqt",1100,850)
c_npe_vs_hqt.cd()
h2_npe_vs_hqt.Draw("colz")
h2_npe_vs_hqt.GetZaxis().SetRangeUser(1,h2_npe_vs_hqt.GetMaximum())
h2_npe_vs_hqt.GetXaxis().SetTitle("log_{10}(NPE)")
h2_npe_vs_hqt.GetYaxis().SetTitle("log_{10}(Homogenized Q_{tot})")
h2_npe_vs_hqt.GetZaxis().SetTitle("Number of Events")
h2_npe_vs_hqt.SetTitle("")
gPad.SetLogz()
gPad.SetRightMargin(0.15)
line.Draw("same")
line.SetLineStyle(9)
c_npe_vs_hqt.SaveAs('npe_vs_hqtot_2d_{:d}events.pdf'.format(int(h2_npe_vs_hqt.GetEntries())))

# h2 zoomed in on highQ region
line2 = ROOT.TLine(3.5,3.5,6,6)
c_npe_vs_hqt_highq = ROOT.TCanvas("c_npe_vs_hqt_highq","c_npe_vs_hqt_highq",1100,850)
c_npe_vs_hqt_highq.cd()
h2_npe_vs_hqt_highq.Draw("colz")
h2_npe_vs_hqt_highq.GetZaxis().SetRangeUser(1,h2_npe_vs_hqt.GetMaximum())
h2_npe_vs_hqt_highq.GetXaxis().SetTitle("log_{10}(NPE)")
h2_npe_vs_hqt_highq.GetYaxis().SetTitle("log_{10}(Homogenized Q_{tot})")
h2_npe_vs_hqt_highq.GetZaxis().SetTitle("Number of Events")
h2_npe_vs_hqt_highq.SetTitle("")
gPad.SetLogz()
gPad.SetRightMargin(0.15)
line2.Draw("same")
line2.SetLineStyle(9)
l_npe_cut = ROOT.TLine(np.log10(25e3),3.5,np.log10(25e3),6)
l_npe_cut.Draw("same")
l_npe_cut.SetLineColor(kRed)
l_npe_cut.SetLineStyle(9)
c_npe_vs_hqt_highq.SaveAs('npe_vs_hqtot_2d_highq_{:d}events.pdf'.format(int(h2_npe_vs_hqt.GetEntries())))

# then the diff plot
c_diff_npe_hqt = ROOT.TCanvas("c_diff_npe_hqt","c_diff_npe_hqt",1100,850)
c_diff_npe_hqt.cd()
h1_diff_npe_hqt.Draw("hist")
h1_diff_npe_hqt.GetXaxis().SetTitle("(NPE-HQtot)/NPE")
h1_diff_npe_hqt.GetYaxis().SetTitle("Number of Events")
h1_diff_npe_hqt.SetTitle("NPE vs HQtot for log_{10}(NPE)>4")
gPad.SetLogy()
c_diff_npe_hqt.SaveAs('npe_vs_hqtot_diff_{:d}events.pdf'.format(int(h2_npe_vs_hqt.GetEntries())))

