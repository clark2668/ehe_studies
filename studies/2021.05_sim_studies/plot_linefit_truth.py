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
cog_var = ['HitStatisticsValues', 'cog_z']
if which_one is 'original':
    charge_var = ['EHEPortiaEventSummarySRT', 'bestNPEbtw']
    zenith_var = ['EHEOpheliaParticleSRT_ImpLF', 'zenith']
    fit_var = ['EHEOpheliaSRT_ImpLF', 'fitQuality']
    speed_var = None
    charge_label = "Portia"
    zenith_label = "Ophelia"
elif which_one is 'new':
    charge_var = ['Homogenized_QTot', 'value']
    zenith_var = ['LineFit_redo', 'zenith']
    speed_var = ['LineFit_redo', 'speed']
    zenith_var_truth = ['PolyplopiaPrimary', 'zenith']
    charge_label = "HQtot"
    zenith_label = "LineFit"

atmo_flux_model = utils_weights.get_flux_model('H3a_SIBYLL23C', 'nugen')
cr_flux_model = utils_weights.get_flux_model('GaisserH3a', 'corsika')

numu_weighter = utils_weights.get_weighter(numu_file, 'nugen', 3000)
nue_weighter = utils_weights.get_weighter(nue_file, 'nugen', 1000)
cor_weighter = utils_weights.get_weighter(cor_file, 'corsika', 1000)

numu_zenith = numu_weighter.get_column(zenith_var[0], zenith_var[1])
numu_zenith_truth = numu_weighter.get_column(zenith_var_truth[0], zenith_var_truth[1])
numu_npe = numu_weighter.get_column(charge_var[0], charge_var[1])
numu_speed = numu_weighter.get_column(speed_var[0], speed_var[1])
numu_inttypes = numu_weighter.get_column('I3MCWeightDict', 'InteractionType')
numu_cog = numu_weighter.get_column(cog_var[0], cog_var[1])

nue_zenith = nue_weighter.get_column(zenith_var[0], zenith_var[1])
nue_zenith_truth = nue_weighter.get_column(zenith_var_truth[0], zenith_var_truth[1])
nue_npe = nue_weighter.get_column(charge_var[0], charge_var[1])
nue_speed = nue_weighter.get_column(speed_var[0], speed_var[1])
nue_cog = nue_weighter.get_column(cog_var[0], cog_var[1])

cor_zenith = cor_weighter.get_column(zenith_var[0], zenith_var[1])
cor_zenith_truth = cor_weighter.get_column(zenith_var_truth[0], zenith_var_truth[1])
cor_npe = cor_weighter.get_column(charge_var[0], charge_var[1])
cor_speed = cor_weighter.get_column(speed_var[0], speed_var[1])
cor_cog = cor_weighter.get_column(cog_var[0], cog_var[1])


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

numu_npe_mask = numu_npe > 25000
nue_npe_mask = nue_npe > 25000
cor_npe_mask = cor_npe > 25000


do_speed_vs_zenith = False
if do_speed_vs_zenith:

    # 2D hist, linefit speed vs true direction

    fig = plt.figure(figsize=(27,7))
    bins = [ np.linspace(-1,1,51), np.linspace(0,.5,51)]
    
    # corsika
    ax = fig.add_subplot(131)
    counts, xedges, yedges, im = ax.hist2d(
            np.cos(cor_zenith_truth[cor_npe_mask]), 
            cor_speed[cor_npe_mask],
            bins=bins,
            weights=cor_weights[cor_npe_mask],
            cmap=cmap,
            norm=colors.LogNorm(),
        )
    # im.set_clim(1E-5, 1E3)
    max_count = np.max(counts)
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Events/Year', fontsize=sizer)
    cbar.ax.tick_params(labelsize=sizer) 
    ax.set_xlabel(r'True $\cos(\theta)$', fontsize=sizer)
    ax.set_ylabel('LineFit Speed', fontsize=sizer)
    ax.tick_params(labelsize=sizer)
    ax.set_title('Corsika (H3a)', fontsize=sizer)

    # numu
    ax2 = fig.add_subplot(132)
    counts, xedges, yedges, im = ax2.hist2d(
            np.cos(numu_zenith_truth[numu_npe_mask]), 
            numu_speed[numu_npe_mask],
            bins=bins,
            weights=numu_atmo_weights[numu_npe_mask],
            cmap=cmap,
            norm=colors.LogNorm(),
        )
    # im.set_clim(1E-5, 1E3)
    cbar2 = plt.colorbar(im, ax=ax2)
    cbar2.set_label('Events/Year', fontsize=sizer)
    cbar2.ax.tick_params(labelsize=sizer) 
    ax2.set_xlabel(r'True $\cos(\theta)$', fontsize=sizer)
    ax2.set_ylabel('LineFit Speed', fontsize=sizer)
    ax2.tick_params(labelsize=sizer)
    ax2.set_title('NuMu (H3a+SIBYLL23C)', fontsize=sizer)

    # nue
    ax3 = fig.add_subplot(133)
    counts, xedges, yedges, im = ax3.hist2d(
            np.cos(nue_zenith_truth[nue_npe_mask]), 
            nue_speed[nue_npe_mask],
            bins=bins,
            weights=nue_atmo_weights[nue_npe_mask],
            cmap=cmap,
            norm=colors.LogNorm(),
        )
    # im.set_clim(1E-5, 1E3)
    cbar3 = plt.colorbar(im, ax=ax3)
    cbar3.set_label('Events/Year', fontsize=sizer)
    cbar3.ax.tick_params(labelsize=sizer) 
    ax3.set_xlabel(r'True $\cos(\theta)$', fontsize=sizer)
    ax3.set_ylabel('LineFit Speed', fontsize=sizer)
    ax3.tick_params(labelsize=sizer)
    ax3.set_title('NuE (H3a+SIBYLL23C)', fontsize=sizer)

    plt.tight_layout()
    fig.savefig('plots/speed_vs_truezenith_{}.png'.format(which_one), edgecolor='none', bbox_inches='tight', dpi=300)
    del fig, ax2, ax3

do_speed_vs_cog = True
if do_speed_vs_cog:

    # 2D hist, linefit speed vs cog
    numu_nc_mask = numu_inttypes > 0

    fig = plt.figure(figsize=(27,7))
    bins = [ np.linspace(-600,600,51), np.linspace(0,.5,51)]
    
    # corsika
    ax = fig.add_subplot(131)
    counts, xedges, yedges, im = ax.hist2d(
            cor_cog[cor_npe_mask],
            cor_speed[cor_npe_mask],
            bins=bins,
            weights=cor_weights[cor_npe_mask],
            cmap=cmap,
            norm=colors.LogNorm(),
        )
    # im.set_clim(1E-5, 1E3)
    max_count = np.max(counts)
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Events/Year', fontsize=sizer)
    cbar.ax.tick_params(labelsize=sizer) 
    ax.set_xlabel(r'Center of Gravity Z', fontsize=sizer)
    ax.set_ylabel('LineFit Speed', fontsize=sizer)
    ax.tick_params(labelsize=sizer)
    ax.set_title('Corsika (H3a)', fontsize=sizer)

    # numu
    ax2 = fig.add_subplot(132)
    counts, xedges, yedges, im = ax2.hist2d(
            numu_cog[numu_npe_mask & numu_nc_mask],
            numu_speed[numu_npe_mask & numu_nc_mask],
            bins=bins,
            weights=numu_atmo_weights[numu_npe_mask & numu_nc_mask],
            cmap=cmap,
            norm=colors.LogNorm(),
        )
    # im.set_clim(1E-5, 1E3)
    cbar2 = plt.colorbar(im, ax=ax2)
    cbar2.set_label('Events/Year', fontsize=sizer)
    cbar2.ax.tick_params(labelsize=sizer) 
    ax2.set_xlabel(r'Center of Gravity Z', fontsize=sizer)
    ax2.set_ylabel('LineFit Speed', fontsize=sizer)
    ax2.tick_params(labelsize=sizer)
    ax2.set_title('NuMu (H3a+SIBYLL23C)', fontsize=sizer)

    # nue
    ax3 = fig.add_subplot(133)
    counts, xedges, yedges, im = ax3.hist2d(
            nue_cog[nue_npe_mask],
            nue_speed[nue_npe_mask],
            bins=bins,
            weights=nue_atmo_weights[nue_npe_mask],
            cmap=cmap,
            norm=colors.LogNorm(),
        )
    # im.set_clim(1E-5, 1E3)
    cbar3 = plt.colorbar(im, ax=ax3)
    cbar3.set_label('Events/Year', fontsize=sizer)
    cbar3.ax.tick_params(labelsize=sizer) 
    ax3.set_xlabel(r'Center of Gravity Z', fontsize=sizer)
    ax3.set_ylabel('LineFit Speed', fontsize=sizer)
    ax3.tick_params(labelsize=sizer)
    ax3.set_title('NuE (H3a+SIBYLL23C)', fontsize=sizer)

    plt.tight_layout()
    fig.savefig('plots/speed_vs_cogz_{}.png'.format(which_one), edgecolor='none', bbox_inches='tight', dpi=300)
    del fig, ax2, ax3
