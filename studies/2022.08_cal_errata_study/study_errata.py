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

# for when we want GZK weights
# gzk_flux = fluxes.EHEFlux("ahlers_gzk")
# gzk_partial = partial(gzk_flux, which_species="nue_sum") 
# the_flux = {
#     'flux': gzk_partial,
#     'label': "GZK_Ahlers"
# }

# for astro weights
def astro_flux(energy):
    # flux of mu @ 100 TeV (per-flavor, per-particle flux)
    # then multiply by two to get the nu+nubar sum
    return 2 * 1.44e-18 / 2 * (energy/1e5)**-2.37

the_flux = {
    'flux': astro_flux,
    'label': "Astro_E2"
}


lims = [1E-5, 1E3]

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

ehe_weights = np.asarray([])
ehe_charge = np.asarray([])
ehe_ndoms = np.asarray([])
ehe_czen = np.asarray([])
ehe_err = np.asarray([])
ehe_sat = np.asarray([])
ehe_edet = np.asarray([])

species = ['nue', 'numu', 'nutau', 'mu', 'tau']
energy_levels = ['high_energy', 'very_high_energy']
species = ['mu']
energy_levels = ['high_energy']
for s in species:
    for level in energy_levels:
        readout_file = cfg_file['juliet'][s][level]['file']
        print(f"Working on juliet {s}, {level}")

        the_f = tables.open_file(readout_file)
        charge = the_f.get_node("/Homogenized_QTot").col('value')
        errata = the_f.get_node("/LenCalErrata").col('value')
        windows = the_f.get_node("/LenSatWindows").col('value')
        energy = the_f.get_node("/I3JulietPrimaryParticle").col('energy')
        ndoms = the_f.get_node(f'/CVMultiplicity').col(f'n_hit_doms')
        # czen = np.cos(the_f.get_node(f'/LineFit').col('zenith'))
        czen = np.cos(the_f.get_node(f'/I3JulietPrimaryParticle').col('zenith'))
        weight_dict, prop_matrix, evts_per_file = weighting.get_juliet_weightdict_and_propmatrix(the_f)
        n_gen = cfg_file['juliet'][s][level]['n_files'] * evts_per_file
        weights = weighting.calc_juliet_weight_from_weight_dict_and_prop_matrix(
            weight_dict=weight_dict, prop_matrix=prop_matrix, 
            flux=the_flux['flux'], n_gen=n_gen, livetime=livetime
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

cor_weights = np.asarray([])
cor_charge = np.asarray([])
cor_ndoms = np.asarray([])
cor_czen = np.asarray([])
cor_err = np.asarray([])
cor_sat = np.asarray([])
cor_edet = np.asarray([])

cor_prim_species = 'fe'
# cor_prim_species = 'p'
cor_sets = [f'21962_{cor_prim_species}', f'22023_{cor_prim_species}']

cr_flux_model = 'GaisserH4a'
cr_flux = weighting.get_flux_model(cr_flux_model, 'corsika')

for c in cor_sets:
    print("Working on corsika {}".format(c))
    cor_file = pd.HDFStore(cfg_file['corsika'][c]['file'])
    cor_weighter = weighting.get_weighter( cor_file, 'corsika',
        cfg_file['corsika'][c]['n_files']
        )
    this_cor_charge = cor_weighter.get_column("Homogenized_QTot", "value")
    this_cor_ndoms = cor_weighter.get_column("CVMultiplicity", "n_hit_doms")
    this_cor_czen = np.cos(cor_weighter.get_column("PolyplopiaPrimary", "zenith"))    
    this_cor_err = cor_weighter.get_column("LenCalErrata", "value")
    this_cor_sat = cor_weighter.get_column("LenSatWindows", "value")
    this_cor_weights = cor_weighter.get_weights(cr_flux) * livetime
    this_cor_edet = cor_weighter.get_column("PolyplopiaPrimary", "energy")

    cor_charge = np.concatenate((cor_charge, this_cor_charge))
    cor_ndoms = np.concatenate((cor_ndoms, this_cor_ndoms))
    cor_czen = np.concatenate((cor_czen, this_cor_czen))
    cor_err = np.concatenate((cor_err, this_cor_err))
    cor_sat = np.concatenate((cor_sat, this_cor_sat))
    cor_weights = np.concatenate((cor_weights, this_cor_weights))
    cor_edet = np.concatenate((cor_edet, this_cor_edet))
    
    cor_file.close()
    

e_slices = {
    'all': [1E2, 1E14],
    # 'vlow': [1E5, 1E6],
    # 'low': [1E6, 1E7],
    # 'mid': [1E7, 1E8],
    # 'high': [1E8, 1E9],
    # 'vhigh': [1E9, 1E10],
    # 'vvhigh': [1E10, 1E11],
}
# e_slices = {}

ehe_q_mask = ehe_charge > 2E4
ehe_ndoms_mask = ehe_ndoms > 100
ehe_q_n_mask = np.logical_and(ehe_q_mask, ehe_ndoms_mask)

cor_q_mask = cor_charge > 2E4
cor_ndoms_mask = cor_ndoms > 100
cor_q_n_mask = np.logical_and(cor_q_mask, cor_ndoms_mask)

ehe_e_slice_errata = {}
ehe_e_slice_saturation = {}

cor_e_slice_errata = {}
cor_e_slice_saturation = {}

for e_slice in e_slices:
    slice_range = e_slices[e_slice]
    print(f"Working on e slice {slice_range}")
    
    ehe_emask = np.logical_and(ehe_edet > slice_range[0], ehe_edet < slice_range[1])
    ehe_mask = np.logical_and(ehe_q_n_mask, ehe_emask)
    
    cor_emask = np.logical_and(cor_edet > slice_range[0], cor_edet < slice_range[1])
    cor_mask = np.logical_and(cor_q_n_mask, cor_emask)
    
    #############################    
    #############################
    ###### 1D Hists of Cal Issues and Saturation for various energy slices
    #############################
    #############################
    
    
    # cal errata juliet slice
    h, _, _ = np.histogram2d(
        ehe_charge[ehe_mask], ehe_err[ehe_mask], 
        weights=ehe_weights[ehe_mask], bins=charge_bins
        )
    h2 = h.sum(axis=0) # collapse against the observable (e.g. # errata)
    del h
    h, _ = np.histogram(y_axis_bin_centers, bins=y_axis, weights=h2,
                           #density = True
                           )
    ehe_e_slice_errata[e_slice] = copy.deepcopy(h)
    del h, h2
    
    # cal errata corsika slice
    h, _, _ = np.histogram2d(
        cor_charge[cor_mask], cor_err[cor_mask],
        weights=cor_weights[cor_mask], bins=charge_bins
    )
    h2 = h.sum(axis=0)
    del h
    h, _ = np.histogram( y_axis_bin_centers, bins=y_axis, weights=h2)
    cor_e_slice_errata[e_slice] = copy.deepcopy(h)
    
    
    # num sat doms juliet slice
    h, _, _ = np.histogram2d(
        ehe_charge[ehe_mask], ehe_sat[ehe_mask], 
        weights=ehe_weights[ehe_mask], bins=energy_bins
        )
    h2 = h.sum(axis=0) # collapse against the observable (e.g. # sat windows)
    del h
    h, _ = np.histogram(y_axis_bin_centers, bins=y_axis, weights=h2,
                           #density = True
                           )
    ehe_e_slice_saturation[e_slice] = copy.deepcopy(h)
    del h, h2
    
    # num sat doms corsika slice
    h, _, _ = np.histogram2d(
        cor_charge[cor_mask], cor_sat[cor_mask],
        weights=cor_weights[cor_mask], bins=energy_bins
    )
    h2 = h.sum(axis=0)
    del h
    h, _ = np.histogram( y_axis_bin_centers, bins=y_axis, weights=h2)
    cor_e_slice_saturation[e_slice] = copy.deepcopy(h)
    
    #############################
    #############################
    ###### 2D Histogram of Cal Issues vs Energy and Charge
    #############################
    #############################
    
    
    # calibration errata juliet
    fig, ((ax1, ax2), (ax3,ax4)) = plt.subplots(2,2, figsize=(12,10))
    a, b, c, im = ax1.hist2d(
        ehe_charge[ehe_mask], ehe_err[ehe_mask], weights=ehe_weights[ehe_mask],
        bins = charge_bins, cmap=my_map, norm=norm
    )
    cbar = plt.colorbar(im, ax=ax1, label='Evts / Year')
    im.set_clim(lims)
    ax1.set_xlabel("Charge / PE")
    ax1.set_ylabel("Num DOMs with Errata")
    ax1.set_xscale('log')

    a, b, c, im2 = ax2.hist2d(
        ehe_edet[ehe_mask], ehe_err[ehe_mask], weights=ehe_weights[ehe_mask],
        bins = energy_bins, cmap=my_map, norm=norm
    )
    cbar2 = plt.colorbar(im2, ax=ax2, label='Evts / Year')
    im2.set_clim(lims)
    ax2.set_xlabel("Energy at Detector / GeV")
    ax2.set_ylabel("Num DOMs with Errata")
    ax2.set_xscale('log')

    # saturation windows
    a, b, c, im3 = ax3.hist2d(
        ehe_charge[ehe_mask], ehe_sat[ehe_mask], weights=ehe_weights[ehe_mask],
        bins = charge_bins, cmap=my_map, norm=norm
    )
    cbar = plt.colorbar(im3, ax=ax3, label='Evts / Year')
    im3.set_clim(lims)
    ax3.set_xlabel("Charge / PE")
    ax3.set_ylabel("Num DOMs with Saturation Windows")
    ax3.set_xscale('log')

    a, b, c, im4 = ax4.hist2d(
        ehe_edet[ehe_mask], ehe_sat[ehe_mask], weights=ehe_weights[ehe_mask],
        bins = energy_bins, cmap=my_map, norm=norm
    )
    cbar2 = plt.colorbar(im4, ax=ax4, label='Evts / Year')
    im4.set_clim(lims)
    ax4.set_xlabel("Energy at Detector / GeV")
    ax4.set_ylabel("Num DOMs with Saturation Windows")
    ax4.set_xscale('log')
    fig.suptitle(f"{species}, {the_flux['label']}", size=20)

    fig.tight_layout()
    fig.savefig(f"cal_issues_vs_charge_juliet_{the_flux['label']}_{slice_range[0]:.1e}_{slice_range[1]:.1e}.png", dpi=300)
    del fig, ax1, ax2, ax3, ax4, cbar2

    #############################    
    #############################
    ###### 2D Histogram of Cal Issues vs cos(zen)
    #############################
    #############################
    

    # first juliet
    fig, (ax1, ax2) = plt.subplots(1,2, figsize=(10,5))
    a, b, c, im1 = ax1.hist2d(
        ehe_czen[ehe_mask], ehe_err[ehe_mask], weights=ehe_weights[ehe_mask],
        bins = saterr_vs_czen_bins, cmap=my_map, norm=norm
    )
    cbar = plt.colorbar(im1, ax=ax1, label='Evts / Year')
    im1.set_clim(lims)
    ax1.set_xlabel(r"cos$\theta$")
    ax1.set_ylabel("Num DOMs with Cal Errata")
    
    a, b, c, im2 = ax2.hist2d(
        ehe_czen[ehe_mask], ehe_sat[ehe_mask], weights=ehe_weights[ehe_mask],
        bins = saterr_vs_czen_bins, cmap=my_map, norm=norm
    )
    cbar = plt.colorbar(im2, ax=ax2, label='Evts / Year')
    im1.set_clim(lims)
    ax2.set_xlabel(r"cos$\theta$")
    ax2.set_ylabel("Num DOMs with Saturation Windows")
    
    fig.suptitle(f"{species}, {the_flux['label']}", size=20)
    fig.tight_layout()
    fig.savefig(f"cal_issues_vs_czen_juliet_{the_flux['label']}_{slice_range[0]:.1e}_{slice_range[1]:.1e}.png", dpi=300)
    del fig, ax1, ax2


    # then corsika
    fig, (ax1, ax2) = plt.subplots(1,2, figsize=(10,5))
    a, b, c, im1 = ax1.hist2d(
        cor_czen[cor_mask], cor_err[cor_mask], weights=cor_weights[cor_mask],
        bins = saterr_vs_czen_bins, cmap=my_map, norm=norm
    )
    cbar = plt.colorbar(im1, ax=ax1, label='Evts / Year')
    im1.set_clim(lims)
    ax1.set_xlabel(r"cos$\theta$")
    ax1.set_ylabel("Num DOMs with Cal Errata")
    
    a, b, c, im2 = ax2.hist2d(
        cor_czen[cor_mask], cor_sat[cor_mask], weights=cor_weights[cor_mask],
        bins = saterr_vs_czen_bins, cmap=my_map, norm=norm
    )
    cbar = plt.colorbar(im2, ax=ax2, label='Evts / Year')
    im1.set_clim(lims)
    ax2.set_xlabel(r"cos$\theta$")
    ax2.set_ylabel("Num DOMs with Saturation Windows")
    
    fig.suptitle(f"{cor_sets}, {cr_flux_model}", size=20)
    fig.tight_layout()
    fig.savefig(f"cal_issues_vs_czen_corsika_{cor_prim_species}_{cr_flux_model}_{slice_range[0]:.1e}_{slice_range[1]:.1e}.png", dpi=300)
    del fig, ax1, ax2

    
#############################    
#############################
###### Save the slices for juliet
#############################
#############################


linestyles = ['-', '--', '-.', ':', '-', '--', '-.', ':']
fig, (ax1, ax2) = plt.subplots(1,2,figsize=(10,5))
for ie, e_slice in enumerate(e_slices):
    slice_range = e_slices[e_slice]
    err_slice = ehe_e_slice_errata[e_slice]
    sat_slice = ehe_e_slice_saturation[e_slice]
    ax1.hist(
        y_axis_bin_centers, bins=y_axis, weights=err_slice,
        histtype='step', linewidth=3, linestyle=linestyles[ie],
        label=f"1E{int(np.log10(slice_range[0]))} < Edet <  1E{int(np.log10(slice_range[1]))}",
    )
    ax2.hist(
        y_axis_bin_centers, bins=y_axis, weights=sat_slice,
        histtype='step', linewidth=3, linestyle=linestyles[ie],
        label=f"1E{int(np.log10(slice_range[0]))} < Edet <  1E{int(np.log10(slice_range[1]))}",
    )
ax1.set_xlabel("Num Cal Errata DOMs")    
ax2.set_xlabel("Num Sat DOMs")

for ax in [ax1, ax2]:
    ax.legend()
    ax.set_ylabel("Num Weighted Events")
    ax.set_yscale('log')
    ax.set_ylim([1E-7, 1E1])

fig.suptitle(f"{species}, {the_flux['label']}", size=20)
fig.tight_layout()
fig.savefig(f"cal_issues_energy_slices_juliet_{the_flux['label']}.png")
del fig, ax1, ax2

#############################    
#############################
###### Save the slices for corsika
#############################
#############################

fig, (ax1, ax2) = plt.subplots(1,2,figsize=(10,5))
for ie, e_slice in enumerate(e_slices):
    slice_range = e_slices[e_slice]
    err_slice = cor_e_slice_errata[e_slice]
    sat_slice = cor_e_slice_saturation[e_slice]
    ax1.hist(
        y_axis_bin_centers, bins=y_axis, weights=err_slice,
        histtype='step', linewidth=3, linestyle=linestyles[ie],
        label=f"1E{int(np.log10(slice_range[0]))} < Edet <  1E{int(np.log10(slice_range[1]))}",
    )
    ax2.hist(
        y_axis_bin_centers, bins=y_axis, weights=sat_slice,
        histtype='step', linewidth=3, linestyle=linestyles[ie],
        label=f"1E{int(np.log10(slice_range[0]))} < Edet <  1E{int(np.log10(slice_range[1]))}",
    )
ax1.set_xlabel("Num Cal Errata DOMs")    
ax2.set_xlabel("Num Sat DOMs")

for ax in [ax1, ax2]:
    ax.legend()
    ax.set_ylabel("Num Weighted Events")
    ax.set_yscale('log')
    ax.set_ylim([1E-7, 1E1])

fig.suptitle(f"{species}, {the_flux['label']}", size=20)
fig.tight_layout()
fig.savefig(f"cal_issues_energy_slices_corsika_{cor_prim_species}_{cr_flux_model}.png")


# # bins = np.logspace(4, 10, 50)
# # nsatdoms_bins= np.linspace(0, 100, 11)
# # nsatdoms_slices = {}

# # for iN in range(len(nsatdoms_bins) -1 ):
# #     slice = [nsatdoms_bins[iN], nsatdoms_bins[iN+1] ]
# #     nsatdoms_slices[slice[0]] = slice
# # # nsatdoms_slices = {}
# # # nsatdoms_slices[0] = [0, 1E3]

# # for iN, this_slice in enumerate(nsatdoms_slices):
# #     slice = nsatdoms_slices[this_slice]
# #     print(f"Working on slice {slice}")
# #     nsatdoms_mask = np.logical_and(windows>= slice[0], windows < slice[1])
# #     total_mask = np.logical_and(q_n_mask, nsatdoms_mask)
# #     fig, ax = plt.subplots(1, 1, figsize=(7,5))
# #     a, b, c, im = ax.hist2d(
# #         energy[total_mask], charge[total_mask], weights=weights[total_mask],
# #         bins = [bins, bins],
# #         cmap=my_map, norm=norm   
# #     )
# #     cbar = plt.colorbar(im, ax=ax, label='Evts / Year, Ahlers GZK')
# #     im.set_clim(lims)
# #     ax.set_xlabel("Energy at Detector / GeV")
# #     ax.set_ylabel("HQTot Charge / PE")
# #     ax.set_xscale('log')
# #     ax.set_yscale('log')
# #     ax.set_aspect('equal')
# #     ax.set_title(f"{species}, {int(slice[0])} < NSatDoms < {int(slice[1])}")

# #     fig.tight_layout()
# #     fig.savefig(f"charge_vs_energy_satdoms_{species}_{int(slice[0])}_{int(slice[1])}.png", dpi=300)
# #     del fig, ax
