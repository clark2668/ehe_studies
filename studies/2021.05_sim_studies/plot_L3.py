import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import argparse

import utils_weights
livetime = 86400 * 365

parser = argparse.ArgumentParser()
parser.add_argument("-numu", type=str,
	dest="numu_file", required=True,
	help="paths to numu file")
parser.add_argument("-nue", type=str, 
	dest="nue_file", required=True,
	help="paths to nue file")
parser.add_argument("-cor", type=str, 
	dest="cor_file", required=True,
	help="paths to corsika file")
args = parser.parse_args()

numu_file = pd.HDFStore(args.numu_file)
nue_file = pd.HDFStore(args.nue_file)
cor_file = pd.HDFStore(args.cor_file)

atmo_flux_model = utils_weights.get_flux_model('H3a_SIBYLL23C', 'nugen')
cr_flux_model = utils_weights.get_flux_model('GaisserH3a', 'corsika')

numu_weighter = utils_weights.get_weighter(numu_file, 'nugen', 1000)
nue_weighter = utils_weights.get_weighter(nue_file, 'nugen', 1000)
cor_weighter = utils_weights.get_weighter(cor_file, 'corsika', 1000)

numu_zenith = numu_weighter.get_column('EHEOpheliaParticleSRT_ImpLF', 'zenith')
numu_chisqured = numu_weighter.get_column('EHEOpheliaSRT_ImpLF', 'fitQuality')
numu_npe = numu_weighter.get_column('EHEPortiaEventSummarySRT', 'bestNPEbtw')
numu_hqtot = numu_weighter.get_column('Homogenized_QTot', 'value')
numu_speed = numu_weighter.get_column('LineFit', 'speed')
numu_atmo_weights = numu_weighter.get_weights(atmo_flux_model)
numu_inttypes = numu_weighter.get_column('I3MCWeightDict', 'InteractionType')
numu_fillratio = numu_weighter.get_column('CascadeFillRatio_L3', 'fillratio_from_mean')
numu_y = numu_weighter.get_column('I3MCWeightDict', 'BjorkenY')
numu_xx = numu_weighter.get_column('PolyplopiaPrimary', 'x')
numu_yy = numu_weighter.get_column('PolyplopiaPrimary', 'y')
numu_zz = numu_weighter.get_column('PolyplopiaPrimary', 'z')

nue_zenith = nue_weighter.get_column('EHEOpheliaParticleSRT_ImpLF', 'zenith')
nue_chisqured = nue_weighter.get_column('EHEOpheliaSRT_ImpLF', 'fitQuality')
nue_npe = nue_weighter.get_column('EHEPortiaEventSummarySRT', 'bestNPEbtw')
nue_hqtot = nue_weighter.get_column('Homogenized_QTot', 'value')
nue_speed = nue_weighter.get_column('LineFit', 'speed')
nue_atmo_weights = nue_weighter.get_weights(atmo_flux_model)
nue_fillratio = nue_weighter.get_column('CascadeFillRatio_L3', 'fillratio_from_mean')

cor_zenith = cor_weighter.get_column('EHEOpheliaParticleSRT_ImpLF', 'zenith')
cor_chisqured = cor_weighter.get_column('EHEOpheliaSRT_ImpLF', 'fitQuality')
cor_npe = cor_weighter.get_column('EHEPortiaEventSummarySRT', 'bestNPEbtw')
cor_hqtot = cor_weighter.get_column('Homogenized_QTot', 'value')
cor_speed = cor_weighter.get_column('LineFit', 'speed')
cor_weights = cor_weighter.get_weights(cr_flux_model)
cor_fillratio = cor_weighter.get_column('CascadeFillRatio_L3', 'fillratio_from_mean')

numu_file.close()
nue_file.close()
cor_file.close()

numu_atmo_weights *= livetime
nue_atmo_weights *= livetime
cor_weights *= livetime

# # do some loop
# for i, fill in enumerate(numu_fillratio):
# 	if np.log10(numu_npe[i]) < 5:
# 		continue
# 	if numu_fillratio[i] < 0.95:
# 		continue
# 	if numu_speed[i] > 0.25:
# 		continue
# 	print('On event {}, with npe {:.2f}, fillR {:.2f}, speed {:.2f}, inelast {:.2f}'.format(i, 
# 		np.log10(numu_npe[i]), numu_fillratio[i], numu_speed[i], numu_y[i]))

cmap=plt.cm.plasma
sizer=15
# numu_nc_mask = numu_inttypes < 2

do_event_dist=False
if do_event_dist:
	fig = plt.figure(figsize=(18,7))
	# bins = [np.linspace(0,500,100), np.linspace(4,8,40)]
	bins = 100
	numu_nc_mask = numu_inttypes < 2

	# top down
	ax = fig.add_subplot(121)
	counts, xedges, yedges, im = ax.hist2d(
			numu_xx, 
			numu_yy, 
			bins=bins,
			cmap=cmap,
			norm=colors.LogNorm(),
			cmin=1
			)
	cbar = plt.colorbar(im, ax=ax)
	cbar.set_label('Events', fontsize=sizer)
	ax.set_ylabel(r'Y', fontsize=sizer)
	ax.set_xlabel(r'X', fontsize=sizer)
	ax.tick_params(labelsize=sizer)
	ax.set_title('Top-Down View')

	# side view
	ax2 = fig.add_subplot(122)
	counts, xedges, yedges, im = ax2.hist2d(
			numu_xx, 
			numu_zz, 
			bins=bins,
			cmap=cmap,
			norm=colors.LogNorm(),
			cmin=1
			)
	# im.set_clim(1E-5, 1E3)
	cbar2 = plt.colorbar(im, ax=ax2)
	cbar2.set_label('Events', fontsize=sizer)
	ax2.set_xlabel(r'X', fontsize=sizer)
	ax2.set_ylabel(r'Z', fontsize=sizer)
	ax2.tick_params(labelsize=sizer)
	ax2.set_title('Side view')

	plt.tight_layout()
	fig.savefig('event_dist.png', edgecolor='none', bbox_inches='tight', dpi=300)

do_inelast_study=False
if do_inelast_study:
	fig = plt.figure(figsize=(18,7))
	# bins = [np.linspace(0,500,100), np.linspace(4,8,40)]
	bins = 100
	numu_nc_mask = numu_inttypes < 2

	# corsika
	ax = fig.add_subplot(121)
	counts, xedges, yedges, im = ax.hist2d(
			numu_y[numu_nc_mask], 
			numu_fillratio[numu_nc_mask], 
			bins=bins,
			weights=numu_atmo_weights[numu_nc_mask],
			cmap=cmap,
			norm=colors.LogNorm(),
			)
	# im.set_clim(1E-5, 1E3)
	cbar = plt.colorbar(im, ax=ax)
	cbar.set_label('Events/Year', fontsize=sizer)
	ax.set_ylabel(r'Fill Ratio', fontsize=sizer)
	ax.set_xlabel(r'Bjorken Y', fontsize=sizer)
	ax.tick_params(labelsize=sizer)
	ax.set_title('NuMu (H3a+SIBYLL23C)')

	# numu
	ax2 = fig.add_subplot(122)
	counts, xedges, yedges, im = ax2.hist2d(
			numu_y[numu_nc_mask], 
			numu_speed[numu_nc_mask], 
			bins=bins,
			weights=numu_atmo_weights[numu_nc_mask],
			cmap=cmap,
			norm=colors.LogNorm(),
			)
	# im.set_clim(1E-5, 1E3)
	cbar2 = plt.colorbar(im, ax=ax2)
	cbar2.set_label('Events/Year', fontsize=sizer)
	ax2.set_xlabel(r'Bjorken Y', fontsize=sizer)
	ax2.set_ylabel(r'LineFit Speed', fontsize=sizer)
	ax2.tick_params(labelsize=sizer)
	ax2.set_title('NuMu (H3a+SIBYLL23C)')

	plt.tight_layout()
	fig.savefig('inelast_study.png', edgecolor='none', bbox_inches='tight', dpi=300)

do_npe_chisqured_cut = False
if do_npe_chisqured_cut:

	numu_npe_mask = numu_npe > 25000
	nue_npe_mask = nue_npe > 25000
	cor_npe_mask = cor_npe > 25000


	# 2D hist, original EHE cut (NPE vs Chi-Square)
	fig = plt.figure(figsize=(27,7))
	bins = [np.linspace(0,500,100), np.linspace(4,8,40)]

	# corsika
	ax = fig.add_subplot(131)
	counts, xedges, yedges, im = ax.hist2d(cor_chisqured[cor_npe_mask], 
			np.log10(cor_npe[cor_npe_mask]), 
			bins=bins,
			weights=cor_weights[cor_npe_mask],
			cmap=cmap,
			norm=colors.LogNorm(),
			)
	im.set_clim(1E-5, 1E3)
	max_count = np.max(counts)
	cbar = plt.colorbar(im, ax=ax)
	cbar.set_label('Events/Year', fontsize=sizer)
	ax.set_ylabel(r'log$_{10}$(Portia NPE)', fontsize=sizer)
	ax.set_xlabel(r'Ophelia Reduced Chi-Square', fontsize=sizer)
	ax.tick_params(labelsize=sizer)
	ax.set_title('Corsika (H3a)')

	# numu
	ax2 = fig.add_subplot(132)
	counts, xedges, yedges, im = ax2.hist2d(numu_chisqured[numu_npe_mask], 
			np.log10(numu_npe[numu_npe_mask]), 
			bins=bins,
			weights=numu_atmo_weights[numu_npe_mask],
			cmap=cmap,
			norm=colors.LogNorm(),
			)
	im.set_clim(1E-5, 1E3)
	cbar2 = plt.colorbar(im, ax=ax2)
	cbar2.set_label('Events/Year', fontsize=sizer)
	ax2.set_ylabel(r'log$_{10}$(Portia NPE)', fontsize=sizer)
	ax2.set_xlabel(r'Ophelia Reduced Chi-Square', fontsize=sizer)
	ax2.tick_params(labelsize=sizer)
	ax2.set_title('NuMu (H3a+SIBYLL23C)')

	# nue
	ax3 = fig.add_subplot(133)
	counts, xedges, yedges, im = ax3.hist2d(nue_chisqured[nue_npe_mask], 
			np.log10(nue_npe[nue_npe_mask]), 
			bins=bins,
			weights=nue_atmo_weights[nue_npe_mask],
			cmap=cmap,
			norm=colors.LogNorm(),
			)
	im.set_clim(1E-5, 1E3)
	cbar3 = plt.colorbar(im, ax=ax3)
	cbar3.set_label('Events/Year', fontsize=sizer)
	ax3.set_ylabel(r'log$_{10}$(Portia NPE)', fontsize=sizer)
	ax3.set_xlabel(r'Ophelia Reduced Chi-Square', fontsize=sizer)
	ax3.tick_params(labelsize=sizer)
	ax3.set_title('NuE (H3a+SIBYLL23C)')

	plt.tight_layout()
	fig.savefig('orig_L3_plots.png', edgecolor='none', bbox_inches='tight', dpi=300)
	del fig, ax

do_fillratio_speed_cut = False
if do_fillratio_speed_cut:

	numu_npe_mask = numu_hqtot > 25000
	numu_nc_mask = numu_inttypes < 2
	nue_npe_mask = nue_hqtot > 25000
	cor_npe_mask = cor_hqtot > 25000

	# 2D hist, potential new EHE cut (Hqtot vs LFspeed)
	fig = plt.figure(figsize=(27,7))
	bins = [np.linspace(0,1,50), np.linspace(0,1,50)]
	# bins = 100

	# corsika
	ax = fig.add_subplot(131)
	counts, xedges, yedges, im = ax.hist2d(
			cor_speed[cor_npe_mask], 
			cor_fillratio[cor_npe_mask],
			bins=bins,
			weights=cor_weights[cor_npe_mask],
			cmap=cmap,
			norm=colors.LogNorm(),
			)
	im.set_clim(1E-5, 1E3)
	max_count = np.max(counts)
	cbar = plt.colorbar(im, ax=ax)
	cbar.set_label('Events/Year', fontsize=sizer)
	ax.set_xlabel(r'LineFit Speed', fontsize=sizer)
	ax.set_ylabel(r'Fill Ratio', fontsize=sizer)
	ax.tick_params(labelsize=sizer)
	ax.set_title('Corsika (H3a)')

	# numu
	# masks = [numu_npe_mask, numu_nc_mask]
	# from functools import reduce
	# numu_total_mask = reduce(np.logical_and, makss)
	ax2 = fig.add_subplot(132)
	counts, xedges, yedges, im = ax2.hist2d(
			numu_speed[numu_npe_mask & numu_nc_mask], 
			numu_fillratio[numu_npe_mask & numu_nc_mask],
			bins=bins,
			weights=numu_atmo_weights[numu_npe_mask & numu_nc_mask],
			cmap=cmap,
			norm=colors.LogNorm(),
			)
	im.set_clim(1E-5, 1E3)
	cbar2 = plt.colorbar(im, ax=ax2)
	cbar2.set_label('Events/Year', fontsize=sizer)
	ax2.set_xlabel(r'LineFit Speed', fontsize=sizer)
	ax2.set_ylabel(r'Fill Ratio', fontsize=sizer)	
	ax2.tick_params(labelsize=sizer)
	ax2.set_title('NuMu (H3a+SIBYLL23C)')

	# nue
	ax3 = fig.add_subplot(133)
	counts, xedges, yedges, im = ax3.hist2d(
			nue_speed[nue_npe_mask], 
			nue_fillratio[nue_npe_mask],
			bins=bins,
			weights=nue_atmo_weights[nue_npe_mask],
			cmap=cmap,
			norm=colors.LogNorm(),
			)
	im.set_clim(1E-5, 1E3)
	cbar3 = plt.colorbar(im, ax=ax3)
	cbar3.set_label('Events/Year', fontsize=sizer)
	ax3.set_xlabel(r'LineFit Speed', fontsize=sizer)
	ax3.set_ylabel(r'Fill Ratio', fontsize=sizer)
	ax3.tick_params(labelsize=sizer)
	ax3.set_title('NuE (H3a+SIBYLL23C)')

	plt.tight_layout()
	fig.savefig('new_L3_plots_2.png', edgecolor='none', bbox_inches='tight', dpi=300)
	del fig, ax

do_hqtot_fillratio_cut = False
if do_hqtot_fillratio_cut:

	numu_hqtot_mask = numu_hqtot > 25000
	numu_nc_mask = numu_inttypes > 0
	# numu_nc_mask = numu_inttypes < 2
	# print(numu_nc_mask)
	nue_hqtot_mask = nue_hqtot > 25000
	cor_hqtot_mask = cor_hqtot > 25000

	# 2D hist, potential new EHE cut (Hqtot vs LFspeed)
	fig = plt.figure(figsize=(27,7))
	bins = [np.linspace(0,1.2,70), np.linspace(4,8,40)]
	# bins = 100

	# corsika
	ax = fig.add_subplot(131)
	counts, xedges, yedges, im = ax.hist2d(
			cor_fillratio[cor_hqtot_mask],
			np.log10(cor_hqtot[cor_hqtot_mask]), 
			bins=bins,
			weights=cor_weights[cor_hqtot_mask],
			cmap=cmap,
			norm=colors.LogNorm(),
			)
	# im.set_clim(1E-5, 1E3)
	max_count = np.max(counts)
	cbar = plt.colorbar(im, ax=ax)
	cbar.set_label('Events/Year', fontsize=sizer)
	ax.set_ylabel(r'log$_{10}$(HQtot)', fontsize=sizer)
	# ax.set_xlabel(r'LineFit Speed', fontsize=sizer)
	ax.set_xlabel(r'Fill Ratio', fontsize=sizer)
	ax.tick_params(labelsize=sizer)
	ax.set_title('Corsika (H3a)')

	# numu
	# masks = [numu_hqtot_mask, numu_nc_mask]
	# from functools import reduce
	# numu_total_mask = reduce(np.logical_and, makss)
	ax2 = fig.add_subplot(132)
	counts, xedges, yedges, im = ax2.hist2d(
			numu_fillratio[numu_hqtot_mask & numu_nc_mask],
			np.log10(numu_hqtot[numu_hqtot_mask & numu_nc_mask]), 
			bins=bins,
			weights=numu_atmo_weights[numu_hqtot_mask & numu_nc_mask],
			cmap=cmap,
			norm=colors.LogNorm(),
			)
	# im.set_clim(1E-5, 1E3)
	cbar2 = plt.colorbar(im, ax=ax2)
	cbar2.set_label('Events/Year', fontsize=sizer)
	ax2.set_ylabel(r'log$_{10}$(HQtot)', fontsize=sizer)
	# ax2.set_xlabel(r'LineFit Speed', fontsize=sizer)
	ax2.set_xlabel(r'Fill Ratio', fontsize=sizer)	
	ax2.tick_params(labelsize=sizer)
	ax2.set_title('NuMu (H3a+SIBYLL23C)')

	# nue
	ax3 = fig.add_subplot(133)
	counts, xedges, yedges, im = ax3.hist2d(
			nue_fillratio[nue_hqtot_mask],
			np.log10(nue_hqtot[nue_hqtot_mask]), 
			bins=bins,
			weights=nue_atmo_weights[nue_hqtot_mask],
			cmap=cmap,
			norm=colors.LogNorm(),
			)
	# im.set_clim(1E-5, 1E3)
	cbar3 = plt.colorbar(im, ax=ax3)
	cbar3.set_label('Events/Year', fontsize=sizer)
	ax3.set_ylabel(r'log$_{10}$(HQtot)', fontsize=sizer)
	ax3.set_xlabel(r'LineFit Speed', fontsize=sizer)
	# ax3.set_xlabel(r'Fill Ratio', fontsize=sizer)
	ax3.tick_params(labelsize=sizer)
	ax3.set_title('NuE (H3a+SIBYLL23C)')

	plt.tight_layout()
	fig.savefig('L3_hqtot_vs_fillratio.png', edgecolor='none', bbox_inches='tight', dpi=300)
	del fig, ax

do_hqtot_speed_cut = True
if do_hqtot_speed_cut:

	numu_hqtot_mask = numu_hqtot > 25000
	numu_nc_mask = numu_inttypes > 0
	numu_nc_mask = numu_inttypes < 2
	# print(numu_nc_mask)
	nue_hqtot_mask = nue_hqtot > 25000
	cor_hqtot_mask = cor_hqtot > 25000

	# one dimesional, to see things quantitatively
	fig = plt.figure(figsize=(5,5))
	bins = np.linspace(0 ,0.5, 150)
	ax = fig.add_subplot(111)
	ax.hist(cor_speed[cor_hqtot_mask], weights=cor_weights[cor_hqtot_mask],
		histtype='step', bins=bins, label='Corsika', 
	)
	ax.hist(numu_speed[numu_hqtot_mask & numu_nc_mask], weights=numu_atmo_weights[numu_hqtot_mask & numu_nc_mask],
		histtype='step', bins=bins, label='Atm NuMu CC'
	)
	ax.hist(numu_speed[numu_hqtot_mask & ~numu_nc_mask], weights=numu_atmo_weights[numu_hqtot_mask & ~numu_nc_mask],
		histtype='step', bins=bins, label='Atm NuMu NC'
	)
	ax.hist(nue_speed[nue_hqtot_mask], weights=nue_atmo_weights[nue_hqtot_mask],
		histtype='step', bins=bins, label='Atm NuE'
	)
	ax.set_ylim([1E-4, 1E4])
	ax.set_ylabel(r'Events/Year')#, fontsize=sizer)
	ax.set_xlabel(r'LineFit Speed')#, fontsize=sizer)
	ax.set_yscale('log')
	ax.legend(loc='upper left')
	plt.tight_layout()
	fig.savefig('L3_speed.png', edgecolor='none', bbox_inches='tight')
	del fig, ax


	# 2D hist, potential new EHE cut (Hqtot vs LFspeed)
	fig = plt.figure(figsize=(27,7))
	bins = [np.linspace(0,0.5,50), np.linspace(4,8,40)]
	# bins = 100

	# corsika
	ax = fig.add_subplot(131)
	counts, xedges, yedges, im = ax.hist2d(
			cor_speed[cor_hqtot_mask], 
			np.log10(cor_hqtot[cor_hqtot_mask]), 
			bins=bins,
			weights=cor_weights[cor_hqtot_mask],
			cmap=cmap,
			norm=colors.LogNorm(),
			)
	im.set_clim(1E-5, 1E3)
	max_count = np.max(counts)
	cbar = plt.colorbar(im, ax=ax)
	cbar.set_label('Events/Year', fontsize=sizer)
	ax.set_ylabel(r'log$_{10}$(HQtot)', fontsize=sizer)
	ax.set_xlabel(r'LineFit Speed', fontsize=sizer)
	ax.tick_params(labelsize=sizer)
	ax.set_title('Corsika (H3a)')

	# numu
	# masks = [numu_hqtot_mask, numu_nc_mask]
	# from functools import reduce
	# numu_total_mask = reduce(np.logical_and, makss)
	ax2 = fig.add_subplot(132)
	counts, xedges, yedges, im = ax2.hist2d(
			numu_speed[numu_hqtot_mask & numu_nc_mask], 
			np.log10(numu_hqtot[numu_hqtot_mask & numu_nc_mask]), 
			bins=bins,
			weights=numu_atmo_weights[numu_hqtot_mask & numu_nc_mask],
			cmap=cmap,
			norm=colors.LogNorm(),
			)
	im.set_clim(1E-5, 1E3)
	cbar2 = plt.colorbar(im, ax=ax2)
	cbar2.set_label('Events/Year', fontsize=sizer)
	ax2.set_ylabel(r'log$_{10}$(HQtot)', fontsize=sizer)
	ax2.set_xlabel(r'LineFit Speed', fontsize=sizer)
	ax2.tick_params(labelsize=sizer)
	ax2.set_title('NuMu (H3a+SIBYLL23C)')

	# nue
	ax3 = fig.add_subplot(133)
	counts, xedges, yedges, im = ax3.hist2d(
			nue_speed[nue_hqtot_mask], 
			np.log10(nue_hqtot[nue_hqtot_mask]), 
			bins=bins,
			weights=nue_atmo_weights[nue_hqtot_mask],
			cmap=cmap,
			norm=colors.LogNorm(),
			)
	im.set_clim(1E-5, 1E3)
	cbar3 = plt.colorbar(im, ax=ax3)
	cbar3.set_label('Events/Year', fontsize=sizer)
	ax3.set_ylabel(r'log$_{10}$(HQtot)', fontsize=sizer)
	ax3.set_xlabel(r'LineFit Speed', fontsize=sizer)
	ax3.tick_params(labelsize=sizer)
	ax3.set_title('NuE (H3a+SIBYLL23C)')

	plt.tight_layout()
	fig.savefig('L3_hqtot_vs_speed.png', edgecolor='none', bbox_inches='tight', dpi=300)
	del fig, ax


