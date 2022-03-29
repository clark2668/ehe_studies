import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import argparse
import utils_plot as up

import utils_weights, utils_ehe
livetime_ic86_i = 2721082.25
livetime_ic86_ii = 2804133.57
livetime_ic86_iii = 2948421.03
livetime = livetime_ic86_i + livetime_ic86_ii + livetime_ic86_iii
# livetime = livetime_ic86_ii

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
parser.add_argument("-data", type=str,
    dest="data_file", required=True,
    help="path to the 10% data file")
args = parser.parse_args()

numu_file = pd.HDFStore(args.numu_file)
nue_file = pd.HDFStore(args.nue_file)
cor_file = pd.HDFStore(args.cor_file)

which_one = 'new'
if which_one is 'original':
    charge_var = ['EHEPortiaEventSummarySRT', 'bestNPEbtw']
    zenith_var = ['EHEOpheliaParticleSRT_ImpLF', 'zenith']
    speed_var = ['EHEOpheliaParticleSRT_ImpLF', 'speed']
    charge_label = "Portia"
    zenith_label = "Ophelia"
    speed_label = "Ophelia Speed"
elif which_one is 'new':
    charge_var = ['Homogenized_QTot', 'value']
    zenith_var = ['LineFit_redo', 'zenith']
    speed_var = ['LineFit_redo', 'speed']
    charge_label = "HQtot"
    zenith_label = "LineFit"
    speed_label = "LineFit Speed"

atmo_model = 'H3a_SIBYLL23C'
cr_model = 'GaisserH4a'
atmo_flux_model = utils_weights.get_flux_model(atmo_model, 'nugen')
cr_flux_model = utils_weights.get_flux_model(cr_model, 'corsika')

numu_weighter = utils_weights.get_weighter(numu_file, 'nugen', 3000)
nue_weighter = utils_weights.get_weighter(nue_file, 'nugen', 1000)
cor_weighter = utils_weights.get_weighter(cor_file, 'corsika', 10000)

charge_masks = {}
log10_charge_cut = 4.4

lims = 1E-5, 1E3

numu_zenith = numu_weighter.get_column(zenith_var[0], zenith_var[1])
numu_coszenith = np.cos(numu_zenith)
numu_speed = numu_weighter.get_column(speed_var[0], speed_var[1])
numu_npe = numu_weighter.get_column(charge_var[0], charge_var[1])
charge_masks['numu'] = np.log10(numu_npe)> log10_charge_cut

nue_zenith = nue_weighter.get_column(zenith_var[0], zenith_var[1])
nue_coszenith = np.cos(nue_zenith)
# nue_speed = nue_weighter.get_column(speed_var[0], speed_var[1])
nue_speed = nue_weighter.get_column('LineFit', 'speed')
nue_npe = nue_weighter.get_column(charge_var[0], charge_var[1])
charge_masks['nue'] = np.log10(nue_npe) > log10_charge_cut

cor_zenith = cor_weighter.get_column(zenith_var[0], zenith_var[1])
cor_coszenith = np.cos(cor_zenith)
cor_speed = cor_weighter.get_column(speed_var[0], speed_var[1])
cor_npe = cor_weighter.get_column(charge_var[0], charge_var[1])
charge_masks['cor'] = np.log10(cor_npe)> log10_charge_cut

numu_atmo_weights = numu_weighter.get_weights(atmo_flux_model)
nue_atmo_weights = nue_weighter.get_weights(atmo_flux_model)
cor_weights = cor_weighter.get_weights(cr_flux_model)

def northern_tracks(energy):
    return 1.44e-18 / 2 * (energy/1e5)**-2.2
numu_astro_weights = numu_weighter.get_weights(northern_tracks)

numu_file.close()
nue_file.close()
cor_file.close()

numu_atmo_weights *= livetime
numu_astro_weights *= livetime
nue_atmo_weights *= livetime
cor_weights *= livetime

data_file = pd.HDFStore(args.data_file)
data_npe = np.asarray(data_file.get(charge_var[0]).get(charge_var[1]))
data_npe = np.log10(data_npe)
charge_masks['data'] = data_npe > log10_charge_cut
data_zenith = np.asarray(data_file.get(zenith_var[0]).get(zenith_var[1]))
data_coszenith = np.cos(data_zenith)
data_speed = np.asarray(data_file.get(speed_var[0]).get(speed_var[1]))

do_L2_plot = True
if do_L2_plot:


    # this shows the 1D distributions of the L2 quantities
    # namely the charge and zenith angle

    #########################
    #########################
    ## Charge Histogram
    #########################
    #########################

    fig, (ax, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [2, 1]})
    bins = np.linspace(4.4, 6, 17)

    # histogram and then plot the data
    data, data_bins = np.histogram(data_npe[charge_masks['data']], bins=bins)
    errs = np.sqrt(data)
    binscenters = np.array([0.5 * (data_bins[i] + data_bins[i+1]) for i in range(len(data_bins)-1)])
    ax.errorbar(binscenters, data, yerr=errs, fmt='ko', label='Burn Sample')

    labels = [
            r'Atm $\nu_{\mu}$' + ', {}'.format(atmo_model), 
            r'Atm $\nu_{e}$'+ ', {}'.format(atmo_model), 
            r'Atm $\mu$' + ', {}'.format(cr_model)
            ]

    # histogram the backgrounds
    sim, sim_bins, patches = ax.hist(
        [
            np.log10(numu_npe)[charge_masks['numu']], 
            # np.log10(numu_npe)[charge_masks['numu']], 
            np.log10(nue_npe)[charge_masks['nue']], 
            np.log10(cor_npe)[charge_masks['cor']]
        ],
        weights=[
            numu_atmo_weights[charge_masks['numu']], 
            # numu_astro_weights[charge_masks['numu']],
            nue_atmo_weights[charge_masks['nue']], 
            cor_weights[charge_masks['cor']]
            ],
        label=labels,
        bins=bins,
        stacked=True
    )
    ax.set_ylabel('Events in {:.2f} days'.format(livetime/60/60/24))
    ax.set_yscale('log')
    ax.set_ylim([1E-1, 1E4])
    ax.legend()

    # ratios

    ratio = data/(sim[0]+sim[1]+sim[2])
    ratio_err = ratio * (errs/data)
    ax2.errorbar(binscenters, ratio, yerr=ratio_err, fmt='ko', label='Data/Sim')
    # ax2.plot(binscenters, ratio, 'ko', label='Sim/Data')
    ax2.set_ylabel('Data/Sim')
    ax2.set_xlabel(r'log$_{10}$(Charge)' + '({}, NPE)'.format(charge_label))
    ax2.plot([4, 6], [1, 1], 'k--') 
    ax2.set_ylim([0.5,1.5])

    plt.tight_layout()
    fig.savefig('L2_charge_hist_{}.png'.format(charge_label))

    del fig, ax, ax2

    #########################
    #########################
    ## Zenith Histogram
    #########################
    #########################


    # now for cos(zen)
    fig, (ax, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [2, 1]})
    bins = np.linspace(-1, 1, 20)

    # histogram and then plot the data
    data, data_bins = np.histogram(data_coszenith[charge_masks['data']], bins=bins)
    errs = np.sqrt(data)
    binscenters = np.array([0.5 * (data_bins[i] + data_bins[i+1]) for i in range(len(data_bins)-1)])
    ax.errorbar(binscenters, data, yerr=errs, fmt='ko', label='Burn Sample')

    # histogram the backgrounds
    sim, sim_bins, patches = ax.hist(
        [
            numu_coszenith[charge_masks['numu']], 
            nue_coszenith[charge_masks['nue']], 
            cor_coszenith[charge_masks['cor']]
        ],
        weights=
            [
                numu_atmo_weights[charge_masks['numu']], 
                nue_atmo_weights[charge_masks['nue']], 
                cor_weights[charge_masks['cor']]
            ],
        label=labels,
        bins=bins,
        stacked=True
    )
    ax.set_ylabel('Events in {:.2f} days'.format(livetime/60/60/24))
    ax.set_yscale('log')
    ax.set_ylim([1E-1, 1E4])
    ax.legend(loc='upper left')

    # ratios
    ratio = data/(sim[0]+sim[1]+sim[2])
    ratio_err = ratio * (errs/data)
    ax2.errorbar(binscenters, ratio, yerr=ratio_err, fmt='ko', label='Data/Sim')
    # ax2.plot(binscenters, data/(sim[0]+sim[1]+sim[2]), 'ko', label='Sim/Data')
    ax2.set_ylabel('Data/Sim')
    ax2.set_xlabel(r'$cos \theta$ ({})'.format(zenith_label))
    ax2.plot([-1.1, 1.1], [1, 1], 'k--') 
    ax2.set_ylim([0.5,1.5])

    plt.tight_layout()
    fig.savefig('L2_coszenith_hist_{}.png'.format(zenith_label))
    del fig, ax, ax2

    #########################
    #########################
    ## LineFit Speed Histogram
    #########################
    #########################

    # now for cos(zen)
    fig, (ax, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [2, 1]})
    bins = np.linspace(0,0.5,50)

    # histogram and then plot the data
    data, data_bins = np.histogram(data_speed[charge_masks['data']], bins=bins)
    errs = np.sqrt(data)
    binscenters = np.array([0.5 * (data_bins[i] + data_bins[i+1]) for i in range(len(data_bins)-1)])
    ax.errorbar(binscenters, data, yerr=errs, fmt='ko', label='Burn Sample')

    # histogram the backgrounds
    sim, sim_bins, patches = ax.hist(
        [
            numu_speed[charge_masks['numu']], 
            nue_speed[charge_masks['nue']], 
            cor_speed[charge_masks['cor']]
        ],
        weights=
            [
                numu_atmo_weights[charge_masks['numu']], 
                nue_atmo_weights[charge_masks['nue']], 
                cor_weights[charge_masks['cor']]
            ],
        label=labels,
        bins=bins,
        stacked=True
    )
    ax.set_ylabel('Events in {:.2f} days'.format(livetime/60/60/24))
    ax.set_yscale('log')
    ax.set_ylim([1E-1, 1E4])
    ax.legend(loc='upper left')

    # ratios
    ratio = data/(sim[0]+sim[1]+sim[2])
    ratio_err = ratio * (errs/data)
    ax2.errorbar(binscenters, ratio, yerr=ratio_err, fmt='ko', label='Data/Sim')
    # ax2.plot(binscenters, data/(sim[0]+sim[1]+sim[2]), 'ko', label='Sim/Data')
    ax2.set_ylabel('Data/Sim')
    ax2.set_xlabel('LineFit Speed')
    ax2.plot([0, 0.5], [1, 1], 'k--') 
    ax2.set_ylim([0,2])

    plt.tight_layout()
    fig.savefig('L2_speed_hist_{}.png'.format(zenith_label))
    del fig, ax, ax2

    #########################
    #########################
    ## 2D HQtot vs Zenith histograms
    #########################
    #########################
    labels = {
        'x': r'$cos \theta$ ({})'.format(zenith_label),
        'y': r'log$_{10}$(Charge)' + "({}, NPE)".format(charge_label)
        }


    # the ax order looks moronic, but the panels make more sense later
    fig, ((ax3, ax, ax2), (ax5, ax6, ax4)) = plt.subplots(2, 3, 
        gridspec_kw={
            'height_ratios': [1, 1],
            'width_ratios': [1, 1, 1]},
        figsize=(15,10)
        )
    bins = [np.linspace(-1,1,20), np.linspace(4.4, 6, 17)]
    my_map = plt.cm.plasma

    # Corsika
    # sim_cor, sim_cor_xedges, sim_cor_yedges, sim_cor_im = ax.hist2d(
    #     cor_coszenith, np.log10(cor_npe), weights=cor_weights,
    #     bins=bins,
    #     cmap=my_map,
    #     norm=colors.LogNorm(),
    #     cmin=1E-3
    # )
    # ax.set_title(r'MC: $\mu$ '+ f'{cr_model}')
    # ax.set_ylabel(r'log$_{10}$(Charge)' + "({}, NPE)".format(charge_label))
    # ax.set_xlabel(r'$cos \theta$ ({})'.format(zenith_label))
    # sim_cor_cbar = plt.colorbar(sim_cor_im, ax=ax)
    # sim_cor_im.set_clim(lims)

    # corsika
    sim_cor, sim_cor_edges, sim_cor_yedges, sim_cor_im = up.make_2D_hist(
        ax, cor_coszenith, np.log10(cor_npe), cor_weights,
        bins, my_map, colors.LogNorm(), lims[0],
        labels['x'], labels['y'],
        r'MC: $\mu$ '+ f'{cr_model}'   
    )
    sim_cor_cbar = plt.colorbar(sim_cor_im, ax=ax)
    sim_cor_im.set_clim(lims)

    # atmospheric neutrinos
    sim_atmo, sim_atmo_xedges, sim_atmo_yedges, sim_atmo_im = up.make_2D_hist(
        ax2, 
        np.concatenate([numu_coszenith, nue_coszenith]), 
        np.concatenate([np.log10(numu_npe), np.log10(nue_npe)]),
        np.concatenate([numu_atmo_weights, nue_atmo_weights]),
        bins, my_map, colors.LogNorm(), lims[0],
        labels['x'], labels['y'],
        r'atm $\nu_{\mu} + \nu_{e}$' + f'{atmo_model}'   
    )
    sim_atmo_cbar = plt.colorbar(sim_atmo_im, ax=ax2)
    sim_atmo_im.set_clim(lims)

    # sum over all MC
    sim, sim_xedges, sim_yedges, sim_im = up.make_2D_hist(
        ax3, 
        np.concatenate([numu_coszenith, nue_coszenith, cor_coszenith]), 
        np.concatenate([np.log10(numu_npe), np.log10(nue_npe), np.log10(cor_npe)]),
        np.concatenate([numu_atmo_weights, nue_atmo_weights, cor_weights]),
        bins, my_map, colors.LogNorm(), lims[0],
        labels['x'], labels['y'],
        r'MC: $\mu$ '+ f'{cr_model}' + r', atm $\nu_{\mu} + \nu_{e}$' + f'{atmo_model}'
    )
    sim_cbar = plt.colorbar(sim_im, ax=ax3)
    sim_im.set_clim(lims)

    # burn sample
    data, data_xedges, data_yedges, data_im = up.make_2D_hist(
        ax5, 
        data_coszenith, 
        data_npe,
        np.ones_like(data_npe),
        bins, my_map, colors.LogNorm(), lims[0],
        labels['x'], labels['y'],
        "Data: Burn Sample, {:.2f} days".format(livetime/60/60/24)
    )
    data_cbar = plt.colorbar(data_im, ax=ax5)
    data_im.set_clim(lims)

    ratio = data/sim
    ratio_im = ax6.pcolorfast(data_xedges, data_yedges, ratio.T,
        cmap=plt.cm.plasma,
        # vmin=-1,
        # vmax=1,
        # norm = up.CenteredNorm(vcenter=1, halfrange=0.5)
        # norm=up.MidpointNormalize(
        #     midpoint=1,
        #     vmin=-1, 
        #     vmax=1,
        # )
    )
    ratio_cbar = plt.colorbar(ratio_im, ax=ax6)
    ax6.set_title('Data/Sim')
    ax6.set_xlabel(labels['x'])
    ax6.set_ylabel(labels['y'])
    ratio_im.set_clim(-1, 2)

    fig.tight_layout()
    fig.savefig('L2_charge_vs_coszen_{}_{}.png'.format(charge_label, zenith_label))

    #########################
    #########################
    ## 2D HQtot vs Speed Histograms
    #########################
    #########################
    labels = {
        'x': r'LineFit Speed ({})'.format(speed_label),
        'y': r'log$_{10}$(Charge)' + "({}, NPE)".format(charge_label)
        }


    # the ax order looks moronic, but the panels make more sense later
    fig, ((ax3, ax, ax2), (ax5, ax6, ax4)) = plt.subplots(2, 3, 
        gridspec_kw={
            'height_ratios': [1, 1],
            'width_ratios': [1, 1, 1]},
        figsize=(15,10)
        )
    bins = [np.linspace(0,0.5,50), np.linspace(4.4, 6, 17)]
    my_map = plt.cm.plasma

    # corsika
    sim_cor, sim_cor_edges, sim_cor_yedges, sim_cor_im = up.make_2D_hist(
        ax, cor_speed, np.log10(cor_npe), cor_weights,
        bins, my_map, colors.LogNorm(), lims[0],
        labels['x'], labels['y'],
        r'MC: $\mu$ '+ f'{cr_model}'   
    )
    sim_cor_cbar = plt.colorbar(sim_cor_im, ax=ax)
    sim_cor_im.set_clim(lims)

    # atmospheric neutrinos
    sim_atmo, sim_atmo_xedges, sim_atmo_yedges, sim_atmo_im = up.make_2D_hist(
        ax2, 
        np.concatenate([numu_speed, nue_speed]), 
        np.concatenate([np.log10(numu_npe), np.log10(nue_npe)]),
        np.concatenate([numu_atmo_weights, nue_atmo_weights]),
        # np.concatenate([nue_speed]), 
        # np.concatenate([np.log10(nue_npe)]),
        # np.concatenate([nue_atmo_weights]),
        bins, my_map, colors.LogNorm(), lims[0],
        labels['x'], labels['y'],
        r'atm $\nu_{\mu} + \nu_{e}$' + f'{atmo_model}'   
    )
    sim_atmo_cbar = plt.colorbar(sim_atmo_im, ax=ax2)
    sim_atmo_im.set_clim(lims)

    # sum over all MC
    sim, sim_xedges, sim_yedges, sim_im = up.make_2D_hist(
        ax3, 
        np.concatenate([numu_speed, nue_speed, cor_speed]), 
        np.concatenate([np.log10(numu_npe), np.log10(nue_npe), np.log10(cor_npe)]),
        np.concatenate([numu_atmo_weights, nue_atmo_weights, cor_weights]),
        bins, my_map, colors.LogNorm(), lims[0],
        labels['x'], labels['y'],
        r'MC: $\mu$ '+ f'{cr_model}' + r', atm $\nu_{\mu} + \nu_{e}$' + f'{atmo_model}'
    )
    sim_cbar = plt.colorbar(sim_im, ax=ax3)
    sim_im.set_clim(lims)

    # burn sample
    data, data_xedges, data_yedges, data_im = up.make_2D_hist(
        ax5, 
        data_speed, 
        data_npe,
        np.ones_like(data_npe),
        bins, my_map, colors.LogNorm(), lims[0],
        labels['x'], labels['y'],
        "Data: Burn Sample, {:.2f} days".format(livetime/60/60/24)
    )
    data_cbar = plt.colorbar(data_im, ax=ax5)
    data_im.set_clim(lims)

    ratio = data/sim
    ratio_im = ax6.pcolorfast(data_xedges, data_yedges, ratio.T,
        cmap=plt.cm.plasma,
    )
    ratio_cbar = plt.colorbar(ratio_im, ax=ax6)
    ax6.set_title('Data/Sim')
    ax6.set_xlabel(labels['x'])
    ax6.set_ylabel(labels['y'])
    ratio_im.set_clim(-1, 2)

    fig.tight_layout()
    fig.savefig('L2_charge_vs_speed_{}_{}.png'.format(charge_label, speed_label))


    del fig, ax, ax2, ax3, ax4, ax5, ax6

    # side study for Shigeru, looking at linefit vs zenith angle for any funny business

    #########################
    #########################
    ## Charge Histogram
    #########################
    #########################

    numu_npe_mask = numu_npe > 25000
    nue_npe_mask = nue_npe > 25000
    cor_npe_mask = cor_npe > 25000
    data_mask = data_npe > np.log10(25000)

    labels = {
        'y': r'LineFit Speed ({})'.format(speed_label),
        'x': r'$\cos(\theta)$' + "({})".format(zenith_label)
        }

    # the ax order looks moronic, but the panels make more sense later
    fig, ((ax3, ax, ax2), (ax5, ax6, ax4)) = plt.subplots(2, 3, 
        gridspec_kw={
            'height_ratios': [1, 1],
            'width_ratios': [1, 1, 1]},
        figsize=(15,10)
        )
    bins = [np.linspace(-1, 1, 21), np.linspace(0, 0.5, 51)]
    my_map = plt.cm.plasma

    # corsika
    sim_cor, sim_cor_edges, sim_cor_yedges, sim_cor_im = up.make_2D_hist(
        ax, cor_coszenith[cor_npe_mask], cor_speed[cor_npe_mask], cor_weights[cor_npe_mask],
        bins, my_map, colors.LogNorm(), lims[0],
        labels['x'], labels['y'],
        r'MC: $\mu$ '+ f'{cr_model}'   
    )
    sim_cor_cbar = plt.colorbar(sim_cor_im, ax=ax)
    sim_cor_im.set_clim(lims)

    # atmospheric neutrinos
    sim_atmo, sim_atmo_xedges, sim_atmo_yedges, sim_atmo_im = up.make_2D_hist(
        ax2, 
        np.concatenate([numu_coszenith[numu_npe_mask], nue_coszenith[nue_npe_mask]]),
        np.concatenate([numu_speed[numu_npe_mask], nue_speed[nue_npe_mask]]), 
        np.concatenate([numu_atmo_weights[numu_npe_mask], nue_atmo_weights[nue_npe_mask]]),
        bins, my_map, colors.LogNorm(), lims[0],
        labels['x'], labels['y'],
        r'atm $\nu_{\mu} + \nu_{e}$' + f'{atmo_model}'   
    )
    sim_atmo_cbar = plt.colorbar(sim_atmo_im, ax=ax2)
    sim_atmo_im.set_clim(lims)

    # sum over all MC
    sim, sim_xedges, sim_yedges, sim_im = up.make_2D_hist(
        ax3, 
        np.concatenate([numu_coszenith[numu_npe_mask], nue_coszenith[nue_npe_mask], cor_coszenith[cor_npe_mask]]),
        np.concatenate([numu_speed[numu_npe_mask], nue_speed[nue_npe_mask], cor_speed[cor_npe_mask]]), 
        np.concatenate([numu_atmo_weights[numu_npe_mask], nue_atmo_weights[nue_npe_mask], cor_weights[cor_npe_mask]]),
        bins, my_map, colors.LogNorm(), lims[0],
        labels['x'], labels['y'],
        r'MC: $\mu$ '+ f'{cr_model}' + r', atm $\nu_{\mu} + \nu_{e}$' + f'{atmo_model}'
    )
    sim_cbar = plt.colorbar(sim_im, ax=ax3)
    sim_im.set_clim(lims)

    # burn sample
    data, data_xedges, data_yedges, data_im = up.make_2D_hist(
        ax5, 
        data_coszenith[data_mask],
        data_speed[data_mask], 
        np.ones_like(data_coszenith[data_mask]),
        bins, my_map, colors.LogNorm(), lims[0],
        labels['x'], labels['y'],
        "Data: Burn Sample, {:.2f} days".format(livetime/60/60/24)
    )
    data_cbar = plt.colorbar(data_im, ax=ax5)
    data_im.set_clim(lims)

    ratio = data/sim
    ratio_im = ax6.pcolorfast(data_xedges, data_yedges, ratio.T,
        cmap=plt.cm.plasma,
    )
    ratio_cbar = plt.colorbar(ratio_im, ax=ax6)
    ax6.set_title('Data/Sim')
    ax6.set_xlabel(labels['x'])
    ax6.set_ylabel(labels['y'])
    ratio_im.set_clim(-1, 2)

    fig.tight_layout()
    fig.savefig('L2_speed_vs_zenith_{}_{}.png'.format(charge_label, speed_label))




data_file.close()