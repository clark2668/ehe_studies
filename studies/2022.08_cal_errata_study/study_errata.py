import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import style
import tables
import pandas as pd
import copy
import yaml

style.use('/home/brian/IceCube/ehe/max_tools/EHE_analysis/eheanalysis/ehe.mplstyle')

from functools import partial
from ehefluxes import fluxes
from eheanalysis import weighting, plotting
gzk_flux = fluxes.EHEFlux("ahlers_gzk")
gzk_partial = partial(gzk_flux, which_species="nue_sum") 
lims = [1E-5, 1E-2]

livetime = 60*60*24*365

num_bins = 50
start_bin_center = 0
stop_bin_center = 500
bin_width = (stop_bin_center - start_bin_center) / num_bins
y_axis = np.arange(start_bin_center-bin_width/2, stop_bin_center + bin_width/2, bin_width)
y_axis_bin_centers = plotting.get_bin_centers(y_axis)

charge_bins = [ np.logspace(3, 8, 16), y_axis ]
energy_bins = [ np.logspace(4, 12, 16), y_axis ]
saterr_vs_czen_bins = [ np.linspace(-1, 1, 20), y_axis ]
my_map = plt.cm.plasma
norm = colors.LogNorm()

cfg_file = '../2022.07_cuts/config.yaml'
cfg_file = yaml.safe_load(open(cfg_file))

the_dict = {
    "mu": "mu_high_energy_merged_998files.hdf5",
    "nue": "nue_high_energy_merged_998files.hdf5",
    "numu": "numu_high_energy_merged_999files.hdf5",
    "nutau": "nutau_high_energy_merged_996files.hdf5",
    "tau": "tau_high_energy_merged_999files.hdf5"
}

ehe_weights = np.asarray([])
ehe_charge = np.asarray([])
ehe_ndoms = np.asarray([])
ehe_czen = np.asarray([])
ehe_err = np.asarray([])
ehe_sat = np.asarray([])
ehe_edet = np.asarray([])

species = ['nue', 'numu', 'nutau', 'mu', 'tau']
for s in species:

    readout_file = the_dict[s]
    print(f"Working on juliet {s}")

    the_f = tables.open_file(readout_file)
    charge = the_f.get_node("/Homogenized_QTot").col('value')
    errata = the_f.get_node("/LenCalErrata").col('value')
    windows = the_f.get_node("/LenSatWindows").col('value')
    energy = the_f.get_node("/I3JulietPrimaryParticle").col('energy')
    ndoms = the_f.get_node(f'/CVMultiplicity').col(f'n_hit_doms')
    # czen = np.cos(the_f.get_node(f'/LineFit').col('zenith'))
    czen = np.cos(the_f.get_node(f'/I3JulietPrimaryParticle').col('zenith'))
    weight_dict, prop_matrix, evts_per_file = weighting.get_juliet_weightdict_and_propmatrix(the_f)
    n_gen = cfg_file['juliet'][s]['high_energy']['n_files'] * evts_per_file
    weights = weighting.calc_juliet_weight_from_weight_dict_and_prop_matrix(
        weight_dict=weight_dict, prop_matrix=prop_matrix, 
        flux=gzk_partial, n_gen=n_gen, livetime=livetime
    )
    
    ehe_weights = np.concatenate((ehe_weights, weights))
    ehe_charge = np.concatenate((ehe_charge, charge))
    ehe_ndoms = np.concatenate((ehe_ndoms, ndoms))
    ehe_czen = np.concatenate((ehe_czen, czen))
    ehe_err = np.concatenate((ehe_err, errata))
    ehe_sat = np.concatenate((ehe_sat, windows))
    ehe_edet = np.concatenate((ehe_edet, energy))
    del charge, errata, windows, energy, ndoms, czen, weight_dict, prop_matrix, evts_per_file, n_gen, weights

    the_f.close()

e_slices = {
    'all': [1E2, 1E14],
    # 'low': [1E6, 1E7],
    # 'mid': [1E7, 1E8],
    # 'high': [1E8, 1E9]
}
# e_slices = {}

ehe_q_mask = ehe_charge > 2E4
ehe_ndoms_mask = ehe_ndoms > 100
ehe_q_n_mask = np.logical_and(ehe_q_mask, ehe_ndoms_mask)

ehe_e_slice_errata = {}
ehe_e_slice_saturation = {}

for e_slice in e_slices:
    slice_range = e_slices[e_slice]
    print(f"Working on e slice {slice_range}")
    
    emask = np.logical_and(ehe_edet > slice_range[0], ehe_edet < slice_range[1])
    mask = np.logical_and(ehe_q_n_mask, emask)
    
    h, _, _ = np.histogram2d(
        ehe_charge[mask], ehe_err[mask], 
        weights=ehe_weights[mask], bins=charge_bins
        )
    h2 = h.sum(axis=0) # collapse against the observable (e.g. # errata)
    del h
    h, _ = np.histogram(y_axis_bin_centers, bins=y_axis, weights=h2,
                           #density = True
                           )
    ehe_e_slice_errata[e_slice] = copy.deepcopy(h)
    del h, h2
    
    h, _, _ = np.histogram2d(
        ehe_charge[mask], ehe_sat[mask], 
        weights=ehe_weights[mask], bins=energy_bins
        )
    h2 = h.sum(axis=0) # collapse against the observable (e.g. # sat windows)
    del h
    h, _ = np.histogram(y_axis_bin_centers, bins=y_axis, weights=h2,
                           #density = True
                           )
    ehe_e_slice_saturation[e_slice] = copy.deepcopy(h)
    del h, h2
    
    # calibration errata
    fig, ((ax1, ax2), (ax3,ax4)) = plt.subplots(2,2, figsize=(12,10))
    a, b, c, im = ax1.hist2d(
        ehe_charge[mask], ehe_err[mask], weights=ehe_weights[mask],
        bins = charge_bins, cmap=my_map, norm=norm
    )
    cbar = plt.colorbar(im, ax=ax1, label='Evts / Year')
    im.set_clim(lims)
    ax1.set_xlabel("Charge / PE")
    ax1.set_ylabel("Num DOMs with Errata")
    ax1.set_xscale('log')

    a, b, c, im2 = ax2.hist2d(
        ehe_edet[mask], ehe_err[mask], weights=ehe_weights[mask],
        bins = energy_bins, cmap=my_map, norm=norm
    )
    cbar2 = plt.colorbar(im2, ax=ax2, label='Evts / Year')
    im2.set_clim(lims)
    ax2.set_xlabel("Energy at Detector / GeV")
    ax2.set_ylabel("Num DOMs with Errata")
    ax2.set_xscale('log')

    # saturation windows
    a, b, c, im3 = ax3.hist2d(
        ehe_charge[mask], ehe_sat[mask], weights=ehe_weights[mask],
        bins = charge_bins, cmap=my_map, norm=norm
    )
    cbar = plt.colorbar(im3, ax=ax3, label='Evts / Year')
    im3.set_clim(lims)
    ax3.set_xlabel("Charge / PE")
    ax3.set_ylabel("Num DOMs with Saturation Windows")
    ax3.set_xscale('log')

    a, b, c, im4 = ax4.hist2d(
        ehe_edet[mask], ehe_sat[mask], weights=ehe_weights[mask],
        bins = energy_bins, cmap=my_map, norm=norm
    )
    cbar2 = plt.colorbar(im4, ax=ax4, label='Evts / Year')
    im4.set_clim(lims)
    ax4.set_xlabel("Energy at Detector / GeV")
    ax4.set_ylabel("Num DOMs with Saturation Windows")
    ax4.set_xscale('log')
    fig.suptitle(f"{species}, Ahlers GZK", size=20)

    fig.tight_layout()
    fig.savefig(f"cal_issues_vs_charge_juliet_{slice_range[0]:.1e}_{slice_range[1]:.1e}.png", dpi=300)
    del fig, ax1, ax2, ax3, ax4, cbar2
    
    fig, (ax1, ax2) = plt.subplots(1,2, figsize=(10,5))
    a, b, c, im1 = ax1.hist2d(
        ehe_czen[mask], ehe_err[mask], weights=ehe_weights[mask],
        bins = saterr_vs_czen_bins, cmap=my_map, norm=norm
    )
    cbar = plt.colorbar(im1, ax=ax1, label='Evts / Year')
    im1.set_clim(lims)
    ax1.set_xlabel(r"cos$\theta$")
    ax1.set_ylabel("Num DOMs with Cal Errata")
    
    a, b, c, im2 = ax2.hist2d(
        ehe_czen[mask], ehe_sat[mask], weights=ehe_weights[mask],
        bins = saterr_vs_czen_bins, cmap=my_map, norm=norm
    )
    cbar = plt.colorbar(im2, ax=ax2, label='Evts / Year')
    im1.set_clim(lims)
    ax2.set_xlabel(r"cos$\theta$")
    ax2.set_ylabel("Num DOMs with Saturation Windows")
    
    fig.suptitle(f"{species}, Ahlers GZK", size=20)
    fig.tight_layout()
    fig.savefig(f"cal_issues_vs_czen_juliet_{slice_range[0]:.1e}_{slice_range[1]:.1e}.png", dpi=300)




# linestyles = ['-', '--', '-.', ':']
# fig, (ax1, ax2) = plt.subplots(1,2,figsize=(10,5))
# for ie, e_slice in enumerate(e_slices):
#     slice_range = e_slices[e_slice]
#     err_slice = e_slice_errata[e_slice]
#     sat_slice = e_slice_saturation[e_slice]
#     ax1.hist(
#         y_axis_bin_centers, bins=y_axis, weights=err_slice,
#         histtype='step', linewidth=3, linestyle=linestyles[ie],
#         label=f"1E{int(np.log10(slice_range[0]))} < Edet <  1E{int(np.log10(slice_range[1]))}",
#     )
#     ax2.hist(
#         y_axis_bin_centers, bins=y_axis, weights=sat_slice,
#         histtype='step', linewidth=3, linestyle=linestyles[ie],
#         label=f"1E{int(np.log10(slice_range[0]))} < Edet <  1E{int(np.log10(slice_range[1]))}",
#     )
# ax1.set_xlabel("Num Cal Errata DOMs")    
# ax2.set_xlabel("Num Sat DOMs")

# for ax in [ax1, ax2]:
#     ax.legend()
#     ax.set_ylabel("Num Weighted Events")

# ax1.set_ylim([0, 0.015])
# ax2.set_ylim([0, 0.04])
# ax1.set_xlim([-20, 500])
# ax2.set_xlim([-20, 200])
# fig.suptitle(f"{species}, Ahlers GZK", size=20)
# fig.tight_layout()
# fig.savefig(f"cal_issues_energy_slices_{species}.png")

# bins = np.logspace(4, 10, 50)
# nsatdoms_bins= np.linspace(0, 100, 11)
# nsatdoms_slices = {}

# for iN in range(len(nsatdoms_bins) -1 ):
#     slice = [nsatdoms_bins[iN], nsatdoms_bins[iN+1] ]
#     nsatdoms_slices[slice[0]] = slice
# # nsatdoms_slices = {}
# # nsatdoms_slices[0] = [0, 1E3]

# for iN, this_slice in enumerate(nsatdoms_slices):
#     slice = nsatdoms_slices[this_slice]
#     print(f"Working on slice {slice}")
#     nsatdoms_mask = np.logical_and(windows>= slice[0], windows < slice[1])
#     total_mask = np.logical_and(q_n_mask, nsatdoms_mask)
#     fig, ax = plt.subplots(1, 1, figsize=(7,5))
#     a, b, c, im = ax.hist2d(
#         energy[total_mask], charge[total_mask], weights=weights[total_mask],
#         bins = [bins, bins],
#         cmap=my_map, norm=norm   
#     )
#     cbar = plt.colorbar(im, ax=ax, label='Evts / Year, Ahlers GZK')
#     im.set_clim(lims)
#     ax.set_xlabel("Energy at Detector / GeV")
#     ax.set_ylabel("HQTot Charge / PE")
#     ax.set_xscale('log')
#     ax.set_yscale('log')
#     ax.set_aspect('equal')
#     ax.set_title(f"{species}, {int(slice[0])} < NSatDoms < {int(slice[1])}")

#     fig.tight_layout()
#     fig.savefig(f"charge_vs_energy_satdoms_{species}_{int(slice[0])}_{int(slice[1])}.png", dpi=300)
#     del fig, ax
