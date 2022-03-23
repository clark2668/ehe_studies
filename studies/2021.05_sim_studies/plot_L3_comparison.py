from shutil import which
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

which_one = 'new'
if which_one is 'original':
    charge_var = ['EHEPortiaEventSummarySRT', 'bestNPEbtw']
    zenith_var = ['EHEOpheliaParticleSRT_ImpLF', 'zenith']
    fit_var = ['EHEOpheliaSRT_ImpLF', 'fitQuality', 'Ophelia']
    speed_var = ['EHEOpheliaParticleSRT_ImpLF', 'speed']
    charge_label = "Portia"
    zenith_label = "Ophelia"
elif which_one is 'new':
    charge_var = ['Homogenized_QTot', 'value']
    zenith_var = ['LineFit_redo', 'zenith']
    speed_var = ['LineFit_redo', 'speed']
    # fit_var = ['LineFit_redoQuality', 'value', 'LineFit_Standard']
    # fit_var = ['LineFit_redoQuality_CutFarAway_120', 'value', 'LineFit_CutFarAway_120m']
    fit_var = ['LineFit_redoQuality_CutFarAway_300', 'value', 'LineFit_CutFarAway_300m']
    charge_label = "HQtot"
    zenith_label = "LineFit"

atmo_flux_model = utils_weights.get_flux_model('H3a_SIBYLL23C', 'nugen')
cr_flux_model = utils_weights.get_flux_model('GaisserH3a', 'corsika')

numu_weighter = utils_weights.get_weighter(numu_file, 'nugen', 3000)
nue_weighter = utils_weights.get_weighter(nue_file, 'nugen', 1000)
cor_weighter = utils_weights.get_weighter(cor_file, 'corsika', 1000)

numu_zenith = numu_weighter.get_column(zenith_var[0], zenith_var[1])
numu_chisqured = numu_weighter.get_column(fit_var[0], fit_var[1])
numu_npe = numu_weighter.get_column(charge_var[0], charge_var[1])
numu_speed = numu_weighter.get_column(speed_var[0], speed_var[1])

nue_zenith = nue_weighter.get_column(zenith_var[0], zenith_var[1])
nue_chisqured = nue_weighter.get_column(fit_var[0], fit_var[1])
nue_npe = nue_weighter.get_column(charge_var[0], charge_var[1])
nue_speed = nue_weighter.get_column(speed_var[0], speed_var[1])

cor_zenith = cor_weighter.get_column(zenith_var[0], zenith_var[1])
cor_chisqured = cor_weighter.get_column(fit_var[0], fit_var[1])
cor_npe = cor_weighter.get_column(charge_var[0], charge_var[1])
cor_speed = cor_weighter.get_column(speed_var[0], speed_var[1])

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

    #################################
    #################################
    # 2D hist, original EHE cut (NPE vs Chi-Square)
    #################################
    #################################
    fig = plt.figure(figsize=(27,7))
    bins = [np.linspace(0,500,100), np.linspace(4.4,6.5,20)]
    
    
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
    fig.savefig('plots/L3_npe_vs_fitqual_{}.png'.format(fit_var[2]), edgecolor='none', bbox_inches='tight', dpi=300)
    del fig, ax2, ax3

    #################################
    #################################
    # 1D plots to try and nail down the CDFs for fit quality
    #################################
    #################################

    fig = plt.figure(figsize=(12,4))
    bins = np.linspace(0 ,500, 51)
    
    ax = fig.add_subplot(131)
    ax.hist(cor_chisqured[cor_npe_mask], weights=cor_weights[cor_npe_mask],
        histtype='step', bins=bins, label='Corsika', 
    )
    ax.hist(numu_chisqured[numu_npe_mask], weights=numu_atmo_weights[numu_npe_mask],
        histtype='step', bins=bins, label='NuMu', 
    )
    ax.hist(nue_chisqured[nue_npe_mask], weights=nue_atmo_weights[nue_npe_mask],
        histtype='step', bins=bins, label='NuE', 
    )
    ax.set_ylim([1E-4, 1E4])
    ax.set_ylabel(r'Events/Year')#, fontsize=sizer)
    ax.set_xlabel(r'LineFit Red Chi Square ({})'.format(zenith_label))#, fontsize=sizer)
    ax.set_yscale('log')

    ax2 = fig.add_subplot(132)
    n_cor, bins_cor, patches_cor = ax2.hist(cor_chisqured[cor_npe_mask], weights=cor_weights[cor_npe_mask],
        histtype='step', bins=bins, label='Corsika (Atm)', density=True
    )
    ax2.hist(numu_chisqured[numu_npe_mask], weights=numu_atmo_weights[numu_npe_mask],
        histtype='step', bins=bins, label='NuMu', density=True,
    )
    n_nue, bins_nue, patches_nue = ax2.hist(nue_chisqured[nue_npe_mask], weights=nue_atmo_weights[nue_npe_mask],
        histtype='step', bins=bins, label='NuE (Atm)', density=True,
    )
    ax2.set_ylabel(r'Normalized Counts (Atm Weighting)')#, fontsize=sizer)
    ax2.set_xlabel(r'LineFit Red Chi Square ({})'.format(zenith_label))#, fontsize=sizer)
    # ax2.set_yscale('log')


    bins_cor = (bins_cor[1:] + bins_cor[:-1])/2
    bins_cor_cut = bins_cor < 100
    integral = np.sum(n_cor[bins_cor_cut]) * (bins[1]-bins[0])
    print("fitqual: the cor integral is {}".format(integral))

    bins_nue = (bins_nue[1:] + bins_nue[:-1])/2
    bins_nue_cut = bins_nue > 100
    integral = np.sum(n_nue[bins_nue_cut]) * (bins[1]- bins[0])
    print("fitqual: the nue integral is {}".format(integral))


    ax3 = fig.add_subplot(133)
    n_cor, bins_cor, patches_cor = ax3.hist(cor_chisqured[cor_npe_mask], weights=cor_weights[cor_npe_mask],
        histtype='step', bins=bins, label='Corsika (Atm)', density=True, cumulative=-1
    )
    ax3.hist(numu_chisqured[numu_npe_mask], weights=numu_atmo_weights[numu_npe_mask],
        histtype='step', bins=bins, label='NuMu', density=True, cumulative=1,
    )
    n_nue, bins_nue, patches_nue = ax3.hist(nue_chisqured[nue_npe_mask], weights=nue_atmo_weights[nue_npe_mask],
        histtype='step', bins=bins, label='NuE (Atm)', density=True, cumulative=1,
    )
    ax3.set_ylabel(r'CDF')#, fontsize=sizer)
    ax3.set_xlabel(r'LineFit Red Chi Square ({})'.format(zenith_label))#, fontsize=sizer)

    ax.legend(loc='upper right')
    plt.tight_layout()
    fig.savefig('plots/L3_fitqual_dist_{}.png'.format(fit_var[2]), edgecolor='none', bbox_inches='tight')
    del fig, ax, ax2, ax3


    print("----------------")


    #################################
    #################################
    # 2D hist, potential new EHE cut (Hqtot vs LFspeed)
    #################################
    #################################

    top_val = 5.2
    bot_val = 4.6
    left_val = 0.25
    right_val = 0.27
    slope = (bot_val - top_val)/(right_val - left_val)
    print(slope)
    def new_cut(speed):
        lognpe_cut = 1e30
        if(speed < left_val):
            lognpe_cut = top_val
        elif(speed < right_val):
            lognpe_cut = top_val + slope*(speed - left_val)
        else:
            lognpe_cut = bot_val
        return lognpe_cut
    vec_function = np.vectorize(new_cut)
    xvals = np.linspace(0,0.6,200)
    yvals = vec_function(xvals)



    fig = plt.figure(figsize=(27,7))
    bins = [np.linspace(0,0.5,50), np.linspace(4.4,6.5,20)]
    # bins = 100

    # corsika
    ax = fig.add_subplot(131)
    counts, xedges, yedges, im = ax.hist2d(
            cor_speed[cor_npe_mask], 
            np.log10(cor_npe[cor_npe_mask]), 
            bins=bins,
            weights=cor_weights[cor_npe_mask],
            cmap=cmap,
            norm=colors.LogNorm(),
            )
    ax.plot(xvals, yvals, linewidth=3)
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
            numu_speed[numu_npe_mask], 
            np.log10(numu_npe[numu_npe_mask]), 
            bins=bins,
            weights=numu_atmo_weights[numu_npe_mask],
            cmap=cmap,
            norm=colors.LogNorm(),
            )
    ax2.plot(xvals, yvals, linewidth=3)
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
            nue_speed[nue_npe_mask], 
            np.log10(nue_npe[nue_npe_mask]), 
            bins=bins,
            weights=nue_atmo_weights[nue_npe_mask],
            cmap=cmap,
            norm=colors.LogNorm(),
            )
    ax3.plot(xvals, yvals, linewidth=3)
    im.set_clim(1E-5, 1E3)
    cbar3 = plt.colorbar(im, ax=ax3)
    cbar3.set_label('Events/Year', fontsize=sizer)
    ax3.set_ylabel(r'log$_{10}$(HQtot)', fontsize=sizer)
    ax3.set_xlabel(r'LineFit Speed', fontsize=sizer)
    ax3.tick_params(labelsize=sizer)
    ax3.set_title('NuE (H3a+SIBYLL23C)')

    plt.tight_layout()
    fig.savefig('plots/L3_npe_vs_speed_{}.png'.format(fit_var[2]), edgecolor='none', bbox_inches='tight')
    del fig, ax

    #################################
    #################################
    # 1D plots to try and nail down the CDFs for speed
    #################################
    #################################

    fig = plt.figure(figsize=(12,4))
    bins = np.linspace(0 ,0.5, 51)
    
    ax = fig.add_subplot(131)
    ax.hist(cor_speed[cor_npe_mask], weights=cor_weights[cor_npe_mask],
        histtype='step', bins=bins, label='Corsika',
    )
    ax.hist(numu_speed[numu_npe_mask], weights=numu_atmo_weights[numu_npe_mask],
        histtype='step', bins=bins, label='NuMu',
    )
    ax.hist(nue_speed[nue_npe_mask], weights=nue_atmo_weights[nue_npe_mask],
        histtype='step', bins=bins, label='NuE',
    )
    ax.set_ylim([1E-4, 1E4])
    ax.set_ylabel(r'Events/Year')#, fontsize=sizer)
    ax.set_xlabel(r'LineFit Speed {}'.format(zenith_label))#, fontsize=sizer)
    ax.set_yscale('log')

    ax2 = fig.add_subplot(132)
    n_cor, bins_cor, patches_cor = ax2.hist(cor_speed[cor_npe_mask], weights=cor_weights[cor_npe_mask],
        histtype='step', bins=bins, label='Corsika (Atm)', density=True
    )
    ax2.hist(numu_speed[numu_npe_mask], weights=numu_atmo_weights[numu_npe_mask],
        histtype='step', bins=bins, label='NuMu', density=True
    )
    n_nue, bins_nue, patches_nue = ax2.hist(nue_speed[nue_npe_mask], weights=nue_atmo_weights[nue_npe_mask],
        histtype='step', bins=bins, label='NuE (Atm)', density=True
    )
    ax2.set_ylabel(r'Normalized Counts (Atm Weighting)')#, fontsize=sizer)
    ax2.set_xlabel(r'LineFit Speed {}'.format(zenith_label))#, fontsize=sizer)

    ax3 = fig.add_subplot(133)
    ax3.hist(cor_speed[cor_npe_mask], weights=cor_weights[cor_npe_mask],
        histtype='step', bins=bins, label='Corsika (Atm)', density=True, cumulative=1
    )
    ax3.hist(numu_speed[numu_npe_mask], weights=numu_atmo_weights[numu_npe_mask],
        histtype='step', bins=bins, label='NuMu', density=True, cumulative=-1
    )
    ax3.hist(nue_speed[nue_npe_mask], weights=nue_atmo_weights[nue_npe_mask],
        histtype='step', bins=bins, label='NuE (Atm)', density=True, cumulative=-1
    )
    ax3.set_ylabel(r'CDF')#, fontsize=sizer)
    ax3.set_xlabel(r'LineFit Speed {}'.format(zenith_label))#, fontsize=sizer)


    bins_cor = (bins_cor[1:] + bins_cor[:-1])/2
    bins_cor_cut = bins_cor > 0.26
    integral = np.sum(n_cor[bins_cor_cut]) * (bins_cor[1]-bins_cor[0])
    print("fitspeed: the cor integral is {}".format(integral))

    bins_nue = (bins_nue[1:] + bins_nue[:-1])/2
    bins_nue_cut = bins_nue < 0.26
    integral = np.sum(n_nue[bins_nue_cut]) * (bins_nue[1]- bins_nue[0])
    print("fitspeed: the nue integral is {}".format(integral))


    ax.legend(loc='upper right')
    plt.tight_layout()
    fig.savefig('plots/L3_speed_dist_{}.png'.format(fit_var[2]), 
        edgecolor='none', bbox_inches='tight', dpi=300)
    del fig, ax


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