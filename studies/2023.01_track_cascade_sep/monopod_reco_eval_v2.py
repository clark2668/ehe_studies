import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import style
import tables
import pandas as pd
import copy
import yaml
import numpy as np


from functools import partial
from ehefluxes import fluxes
from eheanalysis import weighting, plotting, cuts
gzk_flux = fluxes.EHEFlux("cosmogenic_ahlers2010_1E18")
gzk_partial = partial(gzk_flux, which_species="nue_sum") 

style.use('/home/brian/IceCube/ehe/ehe_software/ehe_code/EHE_analysis/eheanalysis/ehe.mplstyle')
# style.use('/data/i3home/baclark/IceCube/ehe/ehe_software/ehe_code/EHE_analysis/eheanalysis/ehe.mplstyle')

livetime = 365 * 24 * 60 * 60
print(livetime)

cfg_file = 'config.yaml'
cfg_file = yaml.safe_load(open(cfg_file))

qcut = 21500
charge_var = 'hqtot'
classifier_var = 'speed'
classifier_bins = np.linspace(0, 0.5, 150)
classifier_cut = 0.27
ndoms_cut = 100

which_reco='monopod'
if which_reco is 'splinempe':
    reco_name='EHE_SplineMPE'
elif which_reco is 'monopod':
    reco_name='EHE_Monopod'
elif which_reco is 'linefit':
    reco_name='EHELineFit'

energy_reco = 'EHE_Monopod'

ehe_weights ={
    "nue": np.asarray([]),
    "mu": np.asarray([])
}
ehe_charge = {
    "nue": np.asarray([]),
    "mu": np.asarray([])
}
ehe_ndoms = {
    "nue": np.asarray([]),
    "mu": np.asarray([])
}
ehe_classifier = {
    "nue": np.asarray([]),
    "mu": np.asarray([])
}
ehe_reco_zen = {
    "nue": np.asarray([]),
    "mu": np.asarray([])
}
ehe_reco_azi = {
    "nue": np.asarray([]),
    "mu": np.asarray([])
}
ehe_monopod_zen = {
    "nue": np.asarray([]),
    "mu": np.asarray([])
}
ehe_monopod_azi = {
    "nue": np.asarray([]),
    "mu": np.asarray([])
}
ehe_true_zen = {
    "nue": np.asarray([]),
    "mu": np.asarray([])
}
ehe_true_azi = {
    "nue": np.asarray([]),
    "mu": np.asarray([])
}
ehe_true_e = {
    "nue": np.asarray([]),
    "mu": np.asarray([])
}
ehe_reco_e = {
    "nue": np.asarray([]),
    "mu": np.asarray([])
}
ehe_contaiment = {
    "nue": np.asarray([]),
    "mu": np.asarray([])
}

juliet_species = ["nue"]
juliet_energy_levels = ["high_energy"]
# juliet_energy_levels = []
events_per_file = {
    "nue_high_energy": 600,
    "nue_very_high_energy": 80,
    "mu_high_energy": 150,
    "mu_very_high_energy": 20
}

which_reco_e = 'emequive'
if which_reco_e is 'depe':
    which_reco_e_var = 'DepE'
    xlabel = 'Deposited Energy'
elif which_reco_e is 'emequive':
    which_reco_e_var = 'EMEquivVisDepE'
    xlabel = 'EM Equiv Energy'
which_reco_e_val = 'value'

which_sample = 'cmc_fine'

containment_selection = 'contained'

#############################
# ehe/cosmogenic flux (juliet)
#############################


for s in juliet_species:
    # if 1==1:
    for l in juliet_energy_levels:

        # print(f"Working on juliet {s} {l}")

        # the_f = tables.open_file(cfg_file['nugen']['21218']['file'])
        the_f = tables.open_file(cfg_file['juliet'][s][l][which_sample]['file'])
        # weight_dict, prop_matrix, evts_per_file = weighting.get_juliet_weightdict_and_propmatrix(the_f)

        # evts_per_file = events_per_file[f"{s}_{l}"] # override to fix L2 issue

        charge = the_f.get_node(f"/{cfg_file['variables'][charge_var]['variable']}").col(f"{cfg_file['variables'][charge_var]['value']}")
        ndoms = the_f.get_node(f"/{cfg_file['variables']['ndoms']['variable']}").col(f"{cfg_file['variables']['ndoms']['value']}")
        classifier = the_f.get_node(f"/{cfg_file['variables'][classifier_var]['variable']}").col(f"{cfg_file['variables'][classifier_var]['value']}")
        reco_zen = the_f.get_node(f"/{reco_name}").col("zenith")
        reco_azi = the_f.get_node(f"/{reco_name}").col("azimuth")
        truth_zen = the_f.get_node("/PolyplopiaPrimary").col("zenith")
        truth_azi = the_f.get_node("/PolyplopiaPrimary").col("azimuth")
        # true_e = the_f.get_node(f"/EMEquivVisDepE").col('value')
        true_e = the_f.get_node(f"/{which_reco_e_var}").col(f"{which_reco_e_val}")
        reco_e = the_f.get_node(f"/{energy_reco}").col("energy")
        containment = the_f.get_node("/EHE_Monopod_Containment").col('value')

        # n_gen = cfg_file['juliet'][s][l]['n_files'] * evts_per_file

        # weights = weighting.calc_juliet_weight_from_weight_dict_and_prop_matrix(
        #     weight_dict=weight_dict, prop_matrix=prop_matrix, 
        #     flux=gzk_partial, n_gen=n_gen, livetime=livetime,
        # )
        weights = np.ones_like(charge)

        ehe_weights[s] = np.concatenate((ehe_weights[s], copy.deepcopy(abs(weights))))
        ehe_charge[s] = np.concatenate((ehe_charge[s], copy.deepcopy(charge)))
        ehe_ndoms[s] = np.concatenate((ehe_ndoms[s], copy.deepcopy(ndoms)))
        ehe_classifier[s] = np.concatenate((ehe_classifier[s], copy.deepcopy(classifier)))
        
        ehe_reco_zen[s] = np.concatenate((ehe_reco_zen[s], copy.deepcopy(reco_zen)))
        ehe_reco_azi[s] = np.concatenate((ehe_reco_azi[s], copy.deepcopy(reco_azi)))
        
        ehe_true_zen[s] = np.concatenate((ehe_true_zen[s], copy.deepcopy(truth_zen)))
        ehe_true_azi[s] = np.concatenate((ehe_true_azi[s], copy.deepcopy(truth_azi)))

        ehe_true_e[s] = np.concatenate((ehe_true_e[s], copy.deepcopy(true_e)))
        ehe_reco_e[s] = np.concatenate((ehe_reco_e[s], copy.deepcopy(reco_e)))
        ehe_contaiment[s] = np.concatenate((ehe_contaiment[s], copy.deepcopy(containment)))

        the_f.close()

ehe_mask = {
    "nue": ehe_charge['nue']>qcut,
    "mu": ehe_charge['mu']>qcut
}

# build masks
for f in ehe_mask.keys():
    qcut = 1000000
    q_mask = ehe_charge[f] > qcut
    ndoms_mask = ehe_ndoms[f] > ndoms_cut
    # track_qual_mask = cuts.track_quality_cut(ehe_classifier[f], ehe_charge[f])
    track_mask = ehe_classifier[f] < 0.27
    total_mask = np.logical_and(q_mask, ndoms_mask)
    # total_mask = np.logical_and(total_mask, track_qual_mask)
    total_mask = np.logical_and(total_mask, track_mask)
    if containment_selection is 'uncontained':
        contained_mask = ehe_contaiment[f] < 1
        total_mask = np.logical_and(total_mask, contained_mask)
    elif containment_selection is 'contained':
        contained_mask = ehe_contaiment[f] > 0
        total_mask = np.logical_and(total_mask, contained_mask)
    ehe_mask[f] = total_mask
    # ehe_mask[f] = ehe_charge[f] > 0


cmap=plt.cm.plasma

make_plots = True
if make_plots:

    clims = [1E-5, 1E-2]
    norm = colors.LogNorm()
    czen_bins = np.linspace(-1, 1, 40)
    azi_bins = np.linspace(0,np.pi*2, 40)
    e_bins = np.logspace(np.log10(1E6),np.log10(1E10),40)

    do_e = True
    if do_e:

        # energycontainment_selection
        fig = plt.figure(figsize=(7,5))
        ax = fig.add_subplot(111)
        vals, xedges, yedges, im = plotting.make_2D_hist(
            ax, 
            ehe_true_e['nue'][ehe_mask['nue']],
            ehe_reco_e['nue'][ehe_mask['nue']],
            bins=e_bins, cmap=plt.cm.viridis,
            weights=ehe_weights['nue'][ehe_mask['nue']],
            norm=norm,
            xlabel=xlabel,
            ylabel='Reco Energy [GeV]',
            title=f"{which_reco}"
        )
        cbar = plt.colorbar(im, ax=ax, label='Evts/Year')
        x, y_med, y_lo, y_hi = plotting.find_contours_2D(
            x_values=ehe_true_e['nue'][ehe_mask['nue']],
            y_values=ehe_reco_e['nue'][ehe_mask['nue']],
            xbins=xedges,
        )
        ax.plot(x, y_med, 'r-', label='Median')
        ax.plot(x, y_lo, 'r-.')
        ax.plot(x, y_hi, 'r-.', label='68% contour')
        ax.legend()
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_title(f"{which_sample}, {which_reco_e}, {containment_selection}")

        # im.set_clim(clims)
        fig.tight_layout()
        fig.savefig(f'./figs/reco_energy_monopod_{which_sample}_{which_reco_e}_{containment_selection}.png')
        del fig, ax, im

        # energy
        fig = plt.figure(figsize=(7,5))
        ax = fig.add_subplot(111)
        x, y_med, y_lo, y_hi = plotting.find_contours_2D(
            x_values=ehe_true_e['nue'][ehe_mask['nue']],
            y_values=(ehe_reco_e['nue'][ehe_mask['nue']] - ehe_true_e['nue'][ehe_mask['nue']])/ehe_true_e['nue'][ehe_mask['nue']],
            xbins=xedges,
        )
        ax.set_title(f"{which_sample}, {which_reco_e}, {containment_selection}")
        ax.fill_between(x, y_med, y_hi, color='C0', alpha=0.2)
        ax.fill_between(x,y_med, y_lo, color='C0', alpha=0.2)
        ax.set_ylabel('(Reco-True)/True')
        ax.set_xlabel(xlabel)
        ax.set_xscale('log')
        ax.set_ylim([-0.05,0.15])
        ax.grid()
        ax.hlines(0.,1E6, 1E10, linestyles='--')
        fig.tight_layout()
        fig.savefig(f'./figs/reco_energy_resolution_monopod_{which_sample}_{which_reco_e}_{containment_selection}.png')
        del fig, ax

    do_dir = False
    if do_dir:

        # zenith
        fig = plt.figure(figsize=(7,5))
        ax = fig.add_subplot(111)
        vals, xedges, yedges, im = plotting.make_2D_hist(
            ax, 
            np.cos(ehe_true_zen['nue'][ehe_mask['nue']]), 
            np.cos(ehe_reco_zen['nue'][ehe_mask['nue']]),
            bins=czen_bins, cmap=plt.cm.viridis,
            weights=ehe_weights['nue'][ehe_mask['nue']],
            norm=norm,
            xlabel='True czen',
            ylabel='Reco czen',
            title=f"{which_reco}"
        )
        cbar = plt.colorbar(im, ax=ax, label='Evts/Year')
        x, y_med, y_lo, y_hi = plotting.find_contours_2D(
            x_values=np.cos(ehe_true_zen['nue'][ehe_mask['nue']]),
            y_values=np.cos(ehe_reco_zen['nue'][ehe_mask['nue']]),
            xbins=xedges,
        )
        ax.plot(x, y_med, 'r-', label='Median')
        ax.plot(x, y_lo, 'r-.')
        ax.plot(x, y_hi, 'r-.', label='68% contour')
        ax.legend()

        im.set_clim(clims)
        fig.tight_layout()
        fig.savefig(f'./figs/reco_czen_{which_reco}.png')
        del fig, ax, im

        # delta zenith
        dczen_bins = np.linspace(-0.25,0.25,40)
        if which_reco is 'splinempe':
            dczen_bins = np.linspace(-2,2,40)
        fig = plt.figure(figsize=(7,5))
        ax = fig.add_subplot(111)
        n, b, p = ax.hist(
            np.cos(ehe_true_zen['nue'][ehe_mask['nue']])-np.cos(ehe_reco_zen['nue'][ehe_mask['nue']]),
            weights=ehe_weights['nue'][ehe_mask['nue']],
            bins=dczen_bins,
            histtype='step',
            lw=3,
        )
        m, r1, r2 = plotting.get_median_quantiles(
            np.cos(ehe_true_zen['nue'][ehe_mask['nue']])-np.cos(ehe_reco_zen['nue'][ehe_mask['nue']]),
            weights=ehe_weights['nue'][ehe_mask['nue']]
        )
        maxval = max(n)
        ax.vlines(m, 0, maxval, 'C1', linestyle='--', label='Median {:.2f}'.format(m))
        ax.vlines(r1, 0, maxval, 'C1', linestyle=':', label='68%: [{:.2f}, {:.2f}]'.format(r1,r2))
        ax.vlines(r2, 0, maxval, 'C1', linestyle=':')
        ax.legend()
        ax.set_xlabel('True-Reco Czen'); ax.set_ylabel('Evts/Yr')
        ax.set_title(f"{which_reco}")
        fig.tight_layout()
        fig.savefig(f'./figs/delta_czen_{which_reco}.png')
        del fig, ax


        # azimuth
        fig = plt.figure(figsize=(7,5))
        ax = fig.add_subplot(111)
        vals, xedges, yedges, im = plotting.make_2D_hist(
            ax, 
            ehe_true_azi['nue'][ehe_mask['nue']], 
            ehe_reco_azi['nue'][ehe_mask['nue']],
            bins=azi_bins, cmap=plt.cm.viridis,
            weights=ehe_weights['nue'][ehe_mask['nue']],
            norm=norm,
            xlabel='True azi',
            ylabel='Reco azi',
            title=f"{which_reco}"
        )
        cbar = plt.colorbar(im, ax=ax, label='Evts/Year')
        x, y_med, y_lo, y_hi = plotting.find_contours_2D(
            x_values=ehe_true_azi['nue'][ehe_mask['nue']],
            y_values=ehe_reco_azi['nue'][ehe_mask['nue']],
            xbins=xedges,
        )
        ax.plot(x, y_med, 'r-', label='Median')
        ax.plot(x, y_lo, 'r-.')
        ax.plot(x, y_hi, 'r-.', label='68% contour')
        ax.legend()

        im.set_clim(clims)
        fig.tight_layout()
        fig.savefig(f'./figs/reco_azi_{which_reco}.png')
        del fig, ax, im

        # delta azimuth
        fig = plt.figure(figsize=(7,5))
        ax = fig.add_subplot(111)
        n, b, p = ax.hist(
            np.rad2deg(ehe_true_azi['nue'][ehe_mask['nue']]-ehe_reco_azi['nue'][ehe_mask['nue']]),
            weights=ehe_weights['nue'][ehe_mask['nue']],
            bins=np.linspace(-10,10,40),
            histtype='step',
            lw=3,
        )
        m, r1, r2 = plotting.get_median_quantiles(
            np.rad2deg(ehe_true_azi['nue'][ehe_mask['nue']]-ehe_reco_azi['nue'][ehe_mask['nue']]),
            weights=ehe_weights['nue'][ehe_mask['nue']]
        )
        maxval = max(n)
        ax.vlines(m, 0, maxval, 'C1', linestyle='--', label='Median {:.2f}'.format(m))
        ax.vlines(r1, 0, maxval, 'C1', linestyle=':', label='68%: [{:.2f}, {:.2f}]'.format(r1,r2))
        ax.vlines(r2, 0, maxval, 'C1', linestyle=':')
        ax.legend()
        ax.set_xlabel('True-Reco Azi [deg]'); ax.set_ylabel('Evts/Yr')
        ax.set_title(f"{which_reco}")
        fig.tight_layout()
        fig.savefig(f'./figs/delta_azi_{which_reco}.png')
        del fig, ax
