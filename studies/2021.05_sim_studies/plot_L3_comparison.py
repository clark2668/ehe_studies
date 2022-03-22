import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import argparse

import utils_weights, utils_ehe
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

which_one = 'original'
if which_one is 'original':
    charge_var = ['EHEPortiaEventSummarySRT', 'bestNPEbtw']
    zenith_var = ['EHEOpheliaParticleSRT_ImpLF', 'zenith']
    fit_var = ['EHEOpheliaSRT_ImpLF', 'fitQuality']
    speed_var = None
    charge_label = "Portia"
    zenith_label = "Ophelia"
elif which_one is 'new':
    charge_var = ['Homogenized_QTot', 'value']
    zenith_var = ['LineFit', 'zenith']
    speed_var = ['LineFit', 'speed']
    # fit_var = ['LineFitQuality', 'value']
    fit_var = ['LineFitQuality_CutFarAway', 'value']
    charge_label = "HQtot"
    zenith_label = "LineFit"

atmo_flux_model = utils_weights.get_flux_model('H3a_SIBYLL23C', 'nugen')
cr_flux_model = utils_weights.get_flux_model('GaisserH3a', 'corsika')

numu_weighter = utils_weights.get_weighter(numu_file, 'nugen', 1000)
nue_weighter = utils_weights.get_weighter(nue_file, 'nugen', 1000)
cor_weighter = utils_weights.get_weighter(cor_file, 'corsika', 1000)

numu_zenith = numu_weighter.get_column(zenith_var[0], zenith_var[1])
numu_chisqured = numu_weighter.get_column(fit_var[0], fit_var[1])
numu_npe = numu_weighter.get_column(charge_var[0], charge_var[1])
# numu_speed = numu_weighter.get_column(speed_var[0], speed_var[1])

nue_zenith = nue_weighter.get_column(zenith_var[0], zenith_var[1])
nue_chisqured = nue_weighter.get_column(fit_var[0], fit_var[1])
nue_npe = nue_weighter.get_column(charge_var[0], charge_var[1])
# nue_speed = nue_weighter.get_column(speed_var[0], speed_var[1])

cor_zenith = cor_weighter.get_column(zenith_var[0], zenith_var[1])
cor_chisqured = cor_weighter.get_column(fit_var[0], fit_var[1])
# cor_chisqured = cor_weighter.get_column('LineFit_redoQuality', 'value')
# cor_chisqured = cor_weighter.get_column('LineFit_redoQuality_CutFarAway', 'value')
cor_npe = cor_weighter.get_column(charge_var[0], charge_var[1])
# cor_speed = cor_weighter.get_column(speed_var[0], speed_var[1])

numu_atmo_weights = numu_weighter.get_weights(atmo_flux_model)
nue_atmo_weights = nue_weighter.get_weights(atmo_flux_model)
cor_weights = cor_weighter.get_weights(cr_flux_model)

numu_file.close()
nue_file.close()
cor_file.close()

numu_atmo_weights *= livetime
nue_atmo_weights *= livetime
cor_weights *= livetime

cmap=plt.cm.viridis
sizer=20

do_npe_chisqured_cut = True
if do_npe_chisqured_cut:

    numu_npe_mask = numu_npe > 25000
    nue_npe_mask = nue_npe > 25000
    cor_npe_mask = cor_npe > 25000


    # 2D hist, original EHE cut (NPE vs Chi-Square)
    fig = plt.figure(figsize=(27,7))
    bins = [np.linspace(0,500,100), np.linspace(4,8,40)]
    
    
    vec_function = np.vectorize(utils_ehe.get_lognpecut_by_fitqual)
    xvals = np.linspace(30,500,499)
    yvals = vec_function(xvals)

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
    cbar.ax.tick_params(labelsize=sizer) 
    ax.plot(xvals, yvals, 'k')
    ax.set_ylabel(r'log10({} Charge)'.format(charge_label), fontsize=sizer)
    ax.set_xlabel(r'{} Reduced Chi-Square'.format(zenith_label), fontsize=sizer)
    ax.tick_params(labelsize=sizer)
    ax.set_title('Corsika (H3a)', fontsize=sizer)

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
    ax2.plot(xvals, yvals, 'k')
    cbar2.set_label('Events/Year', fontsize=sizer)
    cbar2.ax.tick_params(labelsize=sizer) 
    ax2.set_ylabel(r'log10({} Charge)'.format(charge_label), fontsize=sizer)
    ax2.set_xlabel(r'{} Reduced Chi-Square'.format(zenith_label), fontsize=sizer)
    ax2.tick_params(labelsize=sizer)
    ax2.set_title('NuMu (H3a+SIBYLL23C)', fontsize=sizer)

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
    ax3.plot(xvals, yvals, 'k')
    cbar3.set_label('Events/Year', fontsize=sizer)
    cbar3.ax.tick_params(labelsize=sizer) 
    ax3.set_ylabel(r'log10({} Charge)'.format(charge_label), fontsize=sizer)
    ax3.set_xlabel(r'{} Reduced Chi-Square'.format(zenith_label), fontsize=sizer)
    ax3.tick_params(labelsize=sizer)
    ax3.set_title('NuE (H3a+SIBYLL23C)', fontsize=sizer)

    plt.tight_layout()
    fig.savefig('plots/L3_npe_vs_fitqual_{}.png'.format(which_one), edgecolor='none', bbox_inches='tight', dpi=300)
    del fig, ax2, ax3

do_chisqured_speed_cut = False
if do_chisqured_speed_cut:
    numu_npe_mask = numu_npe > 25000
    nue_npe_mask = nue_npe > 25000
    cor_npe_mask = cor_npe > 25000


    fig = plt.figure(figsize=(27,7))
    # bins = [np.linspace(0,500,100), np.linspace(4,8,40)]
    bins = [np.linspace(0,500,100), np.linspace(0, 0.5, 50)]

    # corsika
    ax = fig.add_subplot(131)
    counts, xedges, yedges, im = ax.hist2d(
            cor_chisqured[cor_npe_mask], 
            # np.log10(cor_npe[cor_npe_mask]), 
            cor_speed[cor_npe_mask],
            bins=bins,
            weights=cor_weights[cor_npe_mask],
            cmap=cmap,
            norm=colors.LogNorm(),
            )
    im.set_clim(1E-5, 1E3)
    max_count = np.max(counts)
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Events/Year', fontsize=sizer)
    cbar.ax.tick_params(labelsize=sizer) 
    # ax.plot(xvals, yvals, 'k')
    # ax.set_ylabel(r'log10({} Charge)'.format(charge_label), fontsize=sizer)
    ax.set_ylabel(r'LineFit Speed', fontsize=sizer)
    ax.set_xlabel(r'{} Reduced Chi-Square'.format(zenith_label), fontsize=sizer)
    ax.tick_params(labelsize=sizer)
    ax.set_title('Corsika (H3a)', fontsize=sizer)

    # numu
    ax2 = fig.add_subplot(132)
    counts, xedges, yedges, im = ax2.hist2d(
            numu_chisqured[numu_npe_mask], 
            # np.log10(numu_npe[numu_npe_mask]), 
            numu_speed[numu_npe_mask],
            bins=bins,
            weights=numu_atmo_weights[numu_npe_mask],
            cmap=cmap,
            norm=colors.LogNorm(),
            )
    im.set_clim(1E-5, 1E3)
    cbar2 = plt.colorbar(im, ax=ax2)
    ax2.plot(xvals, yvals, 'k')
    cbar2.set_label('Events/Year', fontsize=sizer)
    cbar2.ax.tick_params(labelsize=sizer) 
    ax2.set_ylabel('LineFit Speed', fontsize=sizer)
    # ax2.set_ylabel(r'log10({} Charge)'.format(charge_label), fontsize=sizer)
    ax2.set_xlabel(r'{} Reduced Chi-Square'.format(zenith_label), fontsize=sizer)
    ax2.tick_params(labelsize=sizer)
    ax2.set_title('NuMu (H3a+SIBYLL23C)', fontsize=sizer)

    # nue
    ax3 = fig.add_subplot(133)
    counts, xedges, yedges, im = ax3.hist2d(
            nue_chisqured[nue_npe_mask], 
            nue_speed[nue_npe_mask],
            # np.log10(nue_npe[nue_npe_mask]), 
            bins=bins,
            weights=nue_atmo_weights[nue_npe_mask],
            cmap=cmap,
            norm=colors.LogNorm(),
            )
    im.set_clim(1E-5, 1E3)
    cbar3 = plt.colorbar(im, ax=ax3)
    ax3.plot(xvals, yvals, 'k')
    cbar3.set_label('Events/Year', fontsize=sizer)
    cbar3.ax.tick_params(labelsize=sizer) 
    # ax3.set_ylabel(r'log10({} Charge)'.format(charge_label), fontsize=sizer)
    ax3.set_ylabel('LineFit Speed', fontsize=sizer)
    ax3.set_xlabel(r'{} Reduced Chi-Square'.format(zenith_label), fontsize=sizer)
    ax3.tick_params(labelsize=sizer)
    ax3.set_title('NuE (H3a+SIBYLL23C)', fontsize=sizer)

    plt.tight_layout()
    fig.savefig('plots/L3_speed_vs_fitqual_{}.png'.format(which_one), edgecolor='none', bbox_inches='tight', dpi=300)
    del fig, ax2, ax3