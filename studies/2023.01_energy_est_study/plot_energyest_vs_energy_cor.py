from re import A
import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import pandas as pd
import tables
from eheanalysis import weighting
import simweights
import yaml

n_gen = 10000 * 150
livetime = 86400 * 365 # 1 year

cfg_file = 'config_umd.yaml'
cfg_file = yaml.safe_load(open(cfg_file))

proton_flux = simweights.FixedFractionFlux(fractions={2212: 1,
                                                      1000260560: 0},
                                           basis=simweights.GaisserH4a())
iron_flux = simweights.FixedFractionFlux(fractions={2212: 0,
                                                    1000260560: 1},
                                           basis=simweights.GaisserH4a())
                                        
which_var = 'charge'
which_var = 'millipedeE'

cor_weighters = []
cor_energy = np.asarray([])
cor_energy_est = np.asarray([])
cor_czen = np.asarray([])
cor_ca_x = np.asarray([])
cor_ca_y = np.asarray([])
cor_ca_z = np.asarray([])

corsika_sets = ['22023']
for c in corsika_sets:
    print(f"Working on corsika {c}")
    cor_file = pd.HDFStore(cfg_file['corsika'][c]['file'])
    cor_weighter = simweights.CorsikaWeighter(
        cor_file, cfg_file['corsika'][c]['n_files']
    )
    cor_weighters.append(cor_weighter)
    this_cor_energy_est = cor_weighter.get_column(
        cfg_file['variables'][which_var]['variable'],
        cfg_file['variables'][which_var]['value'],
    )
    thi_cor_energy = cor_weighter.get_column(
        cfg_file['variables']['corsika_E']['variable'],
        cfg_file['variables']['corsika_E']['value']
    )
    this_cor_zen = np.cos(cor_weighter.get_column(
        cfg_file['variables']['splinempe_zenith']['variable'],
        cfg_file['variables']['splinempe_zenith']['value'],
    ))
    this_cor_ca_x = cor_weighter.get_column( 'ClosestApproach', 'x' )
    this_cor_ca_y = cor_weighter.get_column( 'ClosestApproach', 'y' )
    this_cor_ca_z = cor_weighter.get_column( 'ClosestApproach', 'z' )

    cor_energy_est = np.concatenate((cor_energy_est,this_cor_energy_est))
    cor_energy = np.concatenate((cor_energy, thi_cor_energy))
    cor_czen = np.concatenate((cor_czen,this_cor_zen))
    cor_ca_x = np.concatenate((cor_ca_x,this_cor_ca_x))
    cor_ca_y = np.concatenate((cor_ca_y,this_cor_ca_y))
    cor_ca_z = np.concatenate((cor_ca_z,this_cor_ca_z))

cor_weighter = np.sum(cor_weighters)
cor_proton_weights = cor_weighter.get_weights(proton_flux)*livetime
cor_iron_weights = cor_weighter.get_weights(iron_flux)*livetime

slices = {
    # 0: [-1, -0.86],
    # 1: [-0.86, -0.5],
    # 2: [-0.5, 0],
    # 3: [0, 0.5],
    # 4: [-0.5, 0.86],
    # 5: [0.86, 1],
    6: [-0.2, 0.2]
}

e_plot_bins = np.linspace(2, 9, 30)

######
#### Distribtion of Closest Approach (particularly z)
######


# fig = plt.figure()
# ax = fig.add_subplot(111)
# bins = np.linspace(-1, 1, 20)
# ax.hist(
#     cor_czen[czen_mask], bins=bins, histtype='step',
#     weights=cor_proton_weights[czen_mask]
# )
# ax.set_yscale('log')
# ax.set_xlabel('Closest Approach Depth')
# ax.set_ylabel('Weighted Events')
# fig.tight_layout()
# fig.savefig(f'./plots/czen_hist.png')
# del ax, fig


fig = plt.figure()
ax = fig.add_subplot(111)
bins = np.linspace(-1000, 1000, 50)
for s in slices:
    mask_1 = cor_czen > slices[s][0]
    mask_2 = cor_czen <= slices[s][1]
    czen_mask = np.logical_and(mask_1, mask_2)
    ax.hist(
        cor_ca_z[czen_mask], bins=bins, histtype='step',
        weights=cor_proton_weights[czen_mask],
        label=f"({slices[s][0]}, {slices[s][1]}]"
    )
ax.legend()
ax.set_yscale('log')
ax.set_xlabel('Closest Approach Depth')
ax.set_ylabel('Weighted Events')
fig.tight_layout()
fig.savefig(f'./plots/ca_z_hist.png')
del ax, fig


# ######
# #### Energy Estimator vs True E
# ######

# fig2 = plt.figure()
# ax2 = fig2.add_subplot(111)
# bins = [np.linspace(4,10,60), e_plot_bins]
# my_map = plt.cm.plasma
# my_norm = colors.LogNorm()
# vals, x_edges, y_edges, im = ax2.hist2d(
#     np.log10(cor_energy), np.log10(cor_energy_est), 
#     weights=cor_proton_weights,
#     bins=bins, cmap=my_map, norm=my_norm
# )

# sim_cbar = plt.colorbar(im, ax=ax2, label='unweighted events')
# use_weights = True
# if use_weights:
#     lims = 1E-3, 1E7
#     im.set_clim(lims)
#     sim_cbar.ax.set_ylabel('weighted events', rotation=270)
    
# ax2.set_xlabel('Primary Energy log10(GeV)')
# ax2.set_ylabel(cfg_file['variables'][which_var]['label'])
# fig2.tight_layout()
# fig2.savefig(f'./plots/{which_var}_vs_energy.png')
# del ax2, fig2, im

# ######
# #### Energy Estimator vs cos(zen)
# ######

# # first juliet
# energyest_vs_czen_bins = [ np.linspace(-1, 1, 20), e_plot_bins ]
# fig, ax = plt.subplots(1,1, figsize=(7,5))
# a, b, c, im1 = ax.hist2d(
#     cor_czen, np.log10(cor_energy_est), 
#     weights=cor_proton_weights,
#     bins = energyest_vs_czen_bins, cmap=my_map, norm=my_norm
# )
# cbar = plt.colorbar(im1, ax=ax, label='Evts / Year')
# if use_weights:
#     im1.set_clim(lims)
# ax.set_xlabel(r"SplineMPE cos$\theta$")
# ax.set_ylabel(cfg_file['variables'][which_var]['label'])
# fig.tight_layout()
# fig.savefig(f"./plots/{which_var}_vs_czen.png", dpi=300)
# del fig, ax