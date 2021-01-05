import h5py
import argparse
import numpy as np
import utils.utils_ehe as ehe_utils # original ehe utilities

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



h1_npe = ROOT.TH1D("h1_npe","h1_npe",40,1,8)
h1_nchan = ROOT.TH1D("h1_nchan","h1_nchan",100,0,500)
h1_fitqual = ROOT.TH1D("h1_fitqual","h1_fitqual",100,0,200)

# atL2 = after L2 cuts
h1_npe_atL2 = ROOT.TH1D("h1_npe_atL2","h1_npe_atL2",40,4,8)
h1_zenith_atL2 = ROOT.TH1D("h1_zenith_atL2","h1_zenith_atL2",20,-1,1)
h1_fitqual_atL2 = ROOT.TH1D("h1_fitqual_atL2","h1_fitqual_atL2",100,0,200)
h2_npe_vs_zenith_atL2 = ROOT.TH2D("h2_npe_vs_zenith_atL2","h2_npe_vs_zenith_atL2",20,-1,1,40,4,8)
h2_npe_vs_fitqual_atL2 = ROOT.TH2D("h2_npe_vs_fitqual_atL2","h2_npe_vs_fitqual_atL2",100,0,200,40,4,8)

# atL3 = after L3 cuts
h1_npe_atL3 = ROOT.TH1D("h1_npe_atL3","h1_npe_atL3",40,4,8)
h1_zenith_atL3 = ROOT.TH1D("h1_zenith_atL3","h1_zenith_atL3",20,-1,1)
h1_fitqual_atL3 = ROOT.TH1D("h1_fitqual_atL3","h1_fitqual_atL3",100,0,200)
h2_npe_vs_zenith_atL3 = ROOT.TH2D("h2_npe_vs_zenith_atL3","h2_npe_vs_zenith_atL3",20,-1,1,40,4,8)
h2_npe_vs_fitqual_atL3 = ROOT.TH2D("h2_npe_vs_fitqual_atL3","h2_npe_vs_fitqual_atL3",100,0,200,40,4,8)

# atL4 = after L4 cuts
h1_npe_atL4 = ROOT.TH1D("h1_npe_atL4","h1_npe_atL4",40,4,8)
h1_zenith_atL4 = ROOT.TH1D("h1_zenith_atL4","h1_zenith_atL4",20,-1,1)
h1_fitqual_atL4 = ROOT.TH1D("h1_fitqual_atL4","h1_fitqual_atL4",100,0,200)
h2_npe_vs_zenith_atL4 = ROOT.TH2D("h2_npe_vs_zenith_atL4","h2_npe_vs_zenith_atL4",20,-1,1,40,4,8)
h2_npe_vs_fitqual_atL4 = ROOT.TH2D("h2_npe_vs_fitqual_atL4","h2_npe_vs_fitqual_atL4",100,0,200,40,4,8)


print_L1=True
print_L2=False
print_L3=False
print_L4=False



for file in files:
	print(file)
	file_in = h5py.File(file, "r")
	data = file_in['data']
	portia_npe = data['portia_npe']
	portia_nchan = data['portia_nchan']
	ophelia_zenith = data['ophelia_zenith']
	ophelia_fitqual = data['ophelia_fitqual']
	
	for event, this_portia_npe in enumerate(portia_npe):
		
		h1_npe.Fill(np.log10(portia_npe[event]))
		h1_nchan.Fill(portia_nchan[event])
		h1_fitqual.Fill(ophelia_fitqual[event])

		# # L2 cuts are a cut on NPE, Nchan, and fit quality
		# if(ehe_utils.pass_L2(portia_npe[event], portia_nchan[event], ophelia_fitqual[event])):
			
		# 	h1_npe_atL2.Fill(np.log10(portia_npe[event]))
		# 	h1_zenith_atL2.Fill(np.cos(ophelia_zenith[event]))
		# 	h1_fitqual_atL2.Fill(ophelia_fitqual[event])
		# 	h2_npe_vs_zenith_atL2.Fill(np.cos(ophelia_zenith[event]),np.log10(portia_npe[event]))
		# 	h2_npe_vs_fitqual_atL2.Fill(ophelia_fitqual[event],np.log10(portia_npe[event]))

		# 	# L3 cuts are cut on fit quality and NPE
		# 	if(ehe_utils.pass_L3(portia_npe[event], portia_nchan[event])):
				
		# 		h1_npe_atL3.Fill(np.log10(portia_npe[event]))
		# 		h1_zenith_atL3.Fill(np.cos(ophelia_zenith[event]))
		# 		h1_fitqual_atL3.Fill(ophelia_fitqual[event])
		# 		h2_npe_vs_zenith_atL3.Fill(np.cos(ophelia_zenith[event]),np.log10(portia_npe[event]))
		# 		h2_npe_vs_fitqual_atL3.Fill(ophelia_fitqual[event],np.log10(portia_npe[event]))

		# 		# L4 cuts are cut on cos(zenith) and NPE
		# 		if(ehe_utils.pass_L4(portia_npe[event], ophelia_zenith[event])):
				
		# 			h1_npe_atL4.Fill(np.log10(portia_npe[event]))
		# 			h1_zenith_atL4.Fill(np.cos(ophelia_zenith[event]))
		# 			h1_fitqual_atL4.Fill(ophelia_fitqual[event])
		# 			h2_npe_vs_zenith_atL4.Fill(np.cos(ophelia_zenith[event]),np.log10(portia_npe[event]))
		# 			h2_npe_vs_fitqual_atL4.Fill(ophelia_fitqual[event],np.log10(portia_npe[event]))


	file_in.close()
		

# first, the "all data" plots
if(print_L1):

	# h1 npe all data
	c_npe = ROOT.TCanvas("c_npe","c_npe",1100,850)
	c_npe.cd()
	h1_npe.Draw("hist")
	h1_npe.SetTitle("L1;log_{10}(NPE);Number of Events")
	gPad.SetLogy()
	l_npe_cut = ROOT.TLine(np.log10(25e3),0,np.log10(25e3),h1_npe.GetMaximum())
	l_npe_cut.Draw("same")
	l_npe_cut.SetLineColor(kRed)
	l_npe_cut.SetLineStyle(9)
	c_npe.SaveAs('./ehe_h1_npe_L1_{:d}events.pdf'.format(int(h1_npe.GetEntries())))

	# h1 nchan all data
	c_nchan = ROOT.TCanvas("c_npe","c_npe",1100,850)
	c_nchan.cd()
	h1_nchan.Draw("hist")
	h1_nchan.SetTitle("L1;Nchan;Number of Events")
	gPad.SetLogy()
	l_nchan_cut = ROOT.TLine(100,0,100,h1_nchan.GetMaximum())
	l_nchan_cut.Draw("same")
	l_nchan_cut.SetLineColor(kRed)
	l_nchan_cut.SetLineStyle(9)
	c_nchan.SaveAs('ehe_h1_nchan_L1_{:d}events.pdf'.format(int(h1_nchan.GetEntries())))

	# h1 fit qual all data
	c_fitqual = ROOT.TCanvas("c_fitqual","c_fitqual",1100,850)
	c_fitqual.cd()
	h1_fitqual.Draw("hist")
	h1_fitqual.SetTitle("L1;#chi^{2}/NDF Fit Qual;Number of Events")
	gPad.SetLogy()
	l_fitqual_cut = ROOT.TLine(30,0,30,h1_fitqual.GetMaximum())
	l_fitqual_cut.Draw("same")
	l_fitqual_cut.SetLineColor(kRed)
	l_fitqual_cut.SetLineStyle(9)
	c_fitqual.SaveAs('ehe_h1_fitqual_L1_{:d}events.pdf'.format(int(h1_fitqual.GetEntries())))

if(print_L2):

	# h1 npe at L2
	c_npe_L2 = ROOT.TCanvas("c_npe_L2","c_npe_L2",1100,850)
	c_npe_L2.cd()
	h1_npe_atL2.Draw("hist")
	h1_npe_atL2.SetTitle("L2;log_{10}(NPE);Number of Events")
	gPad.SetLogy()
	c_npe_L2.SaveAs('ehe_h1_npe_L2_{:d}events.pdf'.format(int(h1_npe_atL2.GetEntries())))

	# h1 nchan at L2
	c_zenith_L2 = ROOT.TCanvas("c_zenith_L2","c_zenith_L2",1100,850)
	c_zenith_L2.cd()
	h1_zenith_atL2.Draw("hist")
	h1_zenith_atL2.SetTitle("L2;cos(#theta);Number of Events")
	gPad.SetLogy()
	c_zenith_L2.SaveAs('ehe_h1_zenith_L2_{:d}events.pdf'.format(int(h1_npe_atL2.GetEntries())))

	# h1 fit qual at L2
	c_fitqual_L2 = ROOT.TCanvas("c_fitqual_L2","c_fitqual_L2",1100,850)
	c_fitqual_L2.cd()
	h1_fitqual_atL2.Draw("hist")
	h1_fitqual_atL2.SetTitle("L2;#chi^{2}/NDF Fit Qual;Number of Events")
	gPad.SetLogy()
	c_fitqual_L2.SaveAs('ehe_h1_fitqual_L2_{:d}events.pdf'.format(int(h1_npe_atL2.GetEntries())))

	#h2 npe vs zenith at L2
	c_npe_vs_zenith_L2 = ROOT.TCanvas("c_npe_vs_zenith_L2","c_npe_vs_zenith_L2",1100,850)
	c_npe_vs_zenith_L2.cd()
	h2_npe_vs_zenith_atL2.Draw("colz")
	h2_npe_vs_zenith_atL2.SetTitle("L2;cos(#theta);log_{10}(NPE);Number of Events")
	gPad.SetLogz()
	gPad.SetRightMargin(0.15)
	c_npe_vs_zenith_L2.SaveAs('ehe_h2_npe_vs_zenith_L2_{:d}events.pdf'.format(int(h1_npe_atL2.GetEntries())))

	#h2 npe vs fit qual at L2
	c_npe_vs_fitqual = ROOT.TCanvas("c_npe_vs_fitqual","c_npe_vs_fitqual",1100,850)
	c_npe_vs_fitqual.cd()
	h2_npe_vs_fitqual_atL2.Draw("colz")
	h2_npe_vs_fitqual_atL2.SetTitle("L2;#chi^{2}/NDF Fit Qual;log_{10}(NPE);Number of Events")
	gPad.SetLogz()
	gPad.SetRightMargin(0.15)
	c_npe_vs_fitqual.SaveAs('ehe_h2_npe_vs_fitqual_L2_{:d}events.pdf'.format(int(h1_npe_atL2.GetEntries())))

if(print_L3):

	#first make the L2 plots with the L3 cuts superimposed
	xvals_fitqual = np.linspace(30,200,171)
	func = np.vectorize(ehe_utils.get_lognpecut_by_fitqual)
	yvals_lognpe = func(xvals_fitqual)
	gr_L3_cut = ROOT.TGraph(len(xvals_fitqual),xvals_fitqual,yvals_lognpe)

	#h2 npe vs fit qual at L3
	c_npe_vs_fitqual_L2wcut = ROOT.TCanvas("c_npe_vs_fitqual_L2wcut","c_npe_vs_fitqual_L2wcut",1100,850)
	c_npe_vs_fitqual_L2wcut.cd()
	h2_npe_vs_fitqual_atL2.Draw("colz")
	h2_npe_vs_fitqual_atL2.SetTitle("L2;#chi^{2}/NDF Fit Qual;log_{10}(NPE);Number of Events")
	gr_L3_cut.Draw("sameL")
	gPad.SetLogz()
	gPad.SetRightMargin(0.15)
	c_npe_vs_fitqual_L2wcut.SaveAs('ehe_h2_npe_vs_fitqual_L2wcut_{:d}events.pdf'.format(int(h1_npe_atL2.GetEntries())))

	# h1 npe at L3
	c_npe_L3 = ROOT.TCanvas("c_npe_L3","c_npe_L3",1100,850)
	c_npe_L3.cd()
	h1_npe_atL3.Draw("hist")
	h1_npe_atL3.SetTitle("L3;log_{10}(NPE);Number of Events")
	gPad.SetLogy()
	c_npe_L3.SaveAs('ehe_h1_npe_L3_{:d}events.pdf'.format(int(h1_npe_atL3.GetEntries())))

	# h1 nchan at L3
	c_zenith_L3 = ROOT.TCanvas("c_zenith_L3","c_zenith_L3",1100,850)
	c_zenith_L3.cd()
	h1_zenith_atL3.Draw("hist")
	h1_zenith_atL3.SetTitle("L3;cos(#theta);Number of Events")
	gPad.SetLogy()
	c_zenith_L3.SaveAs('ehe_h1_zenith_L3_{:d}events.pdf'.format(int(h1_npe_atL3.GetEntries())))

	# h1 fit qual at L3
	c_fitqual_L3 = ROOT.TCanvas("c_fitqual_L3","c_fitqual_L3",1100,850)
	c_fitqual_L3.cd()
	h1_fitqual_atL3.Draw("hist")
	h1_fitqual_atL3.SetTitle("L3;#chi^{2}/NDF Fit Qual;Number of Events")
	gPad.SetLogy()
	c_fitqual_L3.SaveAs('ehe_h1_fitqual_L3_{:d}events.pdf'.format(int(h1_npe_atL3.GetEntries())))

	#h2 npe vs zenith at L3
	c_npe_vs_zenith_L3 = ROOT.TCanvas("c_npe_vs_zenith_L3","c_npe_vs_zenith_L3",1100,850)
	c_npe_vs_zenith_L3.cd()
	h2_npe_vs_zenith_atL3.Draw("colz")
	h2_npe_vs_zenith_atL3.SetTitle("L3;cos(#theta);log_{10}(NPE);Number of Events")
	gPad.SetLogz()
	gPad.SetRightMargin(0.15)
	c_npe_vs_zenith_L3.SaveAs('ehe_h2_npe_vs_zenith_L3_{:d}events.pdf'.format(int(h1_npe_atL3.GetEntries())))

	#h2 npe vs fit qual at L3
	c_npe_vs_fitqual_L3 = ROOT.TCanvas("c_npe_vs_fitqual_L3","c_npe_vs_fitqual_L3",1100,850)
	c_npe_vs_fitqual_L3.cd()
	h2_npe_vs_fitqual_atL3.Draw("colz")
	h2_npe_vs_fitqual_atL3.SetTitle("L3;#chi^{2}/NDF Fit Qual;log_{10}(NPE);Number of Events")
	gPad.SetLogz()
	gPad.SetRightMargin(0.15)
	c_npe_vs_fitqual_L3.SaveAs('ehe_h2_npe_vs_fitqual_L3_{:d}events.pdf'.format(int(h1_npe_atL3.GetEntries())))

if(print_L4):

	#first make the L3 plots with the L4 cuts superimposed
	xvals_zenith = np.radians(np.linspace(0,180,181))
	func = np.vectorize(ehe_utils.get_lognpecut_by_zenith)
	yvals_lognpe = func(xvals_zenith)
	gr_L4_cut = ROOT.TGraph(len(xvals_zenith),np.cos(xvals_zenith),yvals_lognpe)

	#h2 npe vs cos(zenith) at L3
	c_npe_vs_zenith_L3wcut = ROOT.TCanvas("c_npe_vs_zenith_L3wcut","c_npe_vs_zenith_L3wcut",1100,850)
	c_npe_vs_zenith_L3wcut.cd()
	h2_npe_vs_zenith_atL3.Draw("colz")
	h2_npe_vs_zenith_atL3.SetTitle("L3;cos(#theta);log_{10}(NPE);Number of Events")
	gr_L4_cut.Draw("sameL")
	gPad.SetLogz()
	gPad.SetRightMargin(0.15)
	c_npe_vs_zenith_L3wcut.SaveAs('ehe_h2_npe_vs_zenith_L3wcut_{:d}events.pdf'.format(int(h1_npe_atL3.GetEntries())))

	# h1 npe at L4
	c_npe_L4 = ROOT.TCanvas("c_npe_L4","c_npe_L4",1100,850)
	c_npe_L4.cd()
	h1_npe_atL4.Draw("hist")
	h1_npe_atL4.SetTitle("L4;log_{10}(NPE);Number of Events")
	gPad.SetLogy()
	c_npe_L4.SaveAs('ehe_h1_npe_L4_{:d}events.pdf'.format(int(h1_npe_atL4.GetEntries())))

	# h1 nchan at L4
	c_zenith_L4 = ROOT.TCanvas("c_zenith_L4","c_zenith_L4",1100,850)
	c_zenith_L4.cd()
	h1_zenith_atL4.Draw("hist")
	h1_zenith_atL4.SetTitle("L4;cos(#theta);Number of Events")
	gPad.SetLogy()
	c_zenith_L4.SaveAs('ehe_h1_zenith_L4_{:d}events.pdf'.format(int(h1_npe_atL4.GetEntries())))

	# h1 fit qual at L4
	c_fitqual_L4 = ROOT.TCanvas("c_fitqual_L4","c_fitqual_L4",1100,850)
	c_fitqual_L4.cd()
	h1_fitqual_atL4.Draw("hist")
	h1_fitqual_atL4.SetTitle("L4;#chi^{2}/NDF Fit Qual;Number of Events")
	gPad.SetLogy()
	c_fitqual_L4.SaveAs('ehe_h1_fitqual_L4_{:d}events.pdf'.format(int(h1_npe_atL4.GetEntries())))

	#h2 npe vs zenith at L4
	c_npe_vs_zenith_L4 = ROOT.TCanvas("c_npe_vs_zenith_L4","c_npe_vs_zenith_L4",1100,850)
	c_npe_vs_zenith_L4.cd()
	h2_npe_vs_zenith_atL4.Draw("colz")
	h2_npe_vs_zenith_atL4.SetTitle("L4;cos(#theta);log_{10}(NPE);Number of Events")
	gPad.SetLogz()
	gPad.SetRightMargin(0.15)
	c_npe_vs_zenith_L4.SaveAs('ehe_h2_npe_vs_zenith_L4_{:d}events.pdf'.format(int(h1_npe_atL4.GetEntries())))

	#h2 npe vs fit qual at L4
	c_npe_vs_fitqual_L4 = ROOT.TCanvas("c_npe_vs_fitqual_L4","c_npe_vs_fitqual_L4",1100,850)
	c_npe_vs_fitqual_L4.cd()
	h2_npe_vs_fitqual_atL4.Draw("colz")
	h2_npe_vs_fitqual_atL4.SetTitle("L4;#chi^{2}/NDF Fit Qual;log_{10}(NPE);Number of Events")
	gPad.SetLogz()
	gPad.SetRightMargin(0.15)
	c_npe_vs_fitqual_L4.SaveAs('ehe_h2_npe_vs_fitqual_L4_{:d}events.pdf'.format(int(h1_npe_atL4.GetEntries())))







