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
import helper

livetime = 86400 * 365 # 1 year

import argparse
parser = argparse.ArgumentParser(description='Do energy estimator')
parser.add_argument(
    "-e", type=str,
    dest='estimator', default='charge',
    help='which estimator' )
args = parser.parse_args()

which_var = args.estimator
# which_var = 'charge'
# which_var = 'millipedeE'
# which_var = 'truncatedE'


cfg_file = 'config_umd.yaml'
cfg_file = yaml.safe_load(open(cfg_file))

proton_flux = simweights.FixedFractionFlux(fractions={2212: 1,
                                                      1000260560: 0},
                                           basis=simweights.GaisserH4a())
iron_flux = simweights.FixedFractionFlux(fractions={2212: 0,
                                                    1000260560: 1},
                                           basis=simweights.GaisserH4a())

my_map = plt.cm.plasma
my_norm = colors.LogNorm()
lims = 1E-3, 1E6


cor_weighters = []
cor_energy = np.asarray([])
cor_energy_est = np.asarray([])
cor_czen = np.asarray([])
cor_ca_x = np.asarray([])
cor_ca_y = np.asarray([])
cor_ca_z = np.asarray([])
cor_rlogl = np.asarray([])
cor_speed = np.asarray([])

corsika_sets = ['22023', '21962']
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
    this_cor_ca_x = cor_weighter.get_column( 'CVStatistics', 'cog_x' )
    this_cor_ca_y = cor_weighter.get_column( 'CVStatistics', 'cog_y' )
    this_cor_ca_z = cor_weighter.get_column( 'CVStatistics', 'cog_z' )
    this_cor_rlogl = cor_weighter.get_column('EHE_SplineMPEFitParams', 'rlogl')
    this_cor_speed = cor_weighter.get_column('EHELineFit', 'speed')

    cor_energy_est = np.concatenate((cor_energy_est,this_cor_energy_est))
    cor_energy = np.concatenate((cor_energy, thi_cor_energy))
    cor_czen = np.concatenate((cor_czen,this_cor_zen))
    cor_ca_x = np.concatenate((cor_ca_x,this_cor_ca_x))
    cor_ca_y = np.concatenate((cor_ca_y,this_cor_ca_y))
    cor_ca_z = np.concatenate((cor_ca_z,this_cor_ca_z))
    cor_rlogl = np.concatenate((cor_rlogl, this_cor_rlogl))
    cor_speed = np.concatenate((cor_speed, this_cor_speed))

cor_weighter = np.sum(cor_weighters)
cor_proton_weights = cor_weighter.get_weights(proton_flux)*livetime
cor_iron_weights = cor_weighter.get_weights(iron_flux)*livetime
which_weights = cor_iron_weights
print("Got here")
slices = {
    # 0: [-1, -0.86],
    # 1: [-0.86, -0.5],
    # 2: [-0.5, 0],
    # 3: [0, 0.5],
    # 4: [-0.5, 0.86],
    # 5: [0.86, 1],
    6: [-0.3, 0.3]
}

slices = {}
starter = -1
for i in range(0,10):
    slices[i] = [-1+(i*0.2), -1+((i+1)*0.2)]

e_plot_bins = np.linspace(2, 9, 60)

######
#### Distribtion of Closest Approach (particularly z)
######

fig, axarr =  plt.subplots(5,2, figsize=(8,20))
axarr = axarr.flatten()
bins = np.linspace(-800, 800, 50)
for i, s in enumerate(slices.keys()):
    mask_1 = cor_czen > slices[s][0]
    mask_2 = cor_czen <= slices[s][1]
    czen_mask = np.logical_and(mask_1, mask_2)
    axarr[i].hist(
        cor_ca_z[czen_mask], bins=bins, histtype='step',
        weights=which_weights[czen_mask],
        label=f"({slices[s][0]:.2f}, {slices[s][1]:.2f}]"
    )
    axarr[i].legend()
    axarr[i].set_yscale('log')
    axarr[i].set_xlabel('COG Depth')
    axarr[i].set_ylabel('Weighted Events')
    axarr[i].set_ylim([10,1e7])
fig.tight_layout()
fig.savefig(f'./plots/cog_z_hist.png', dpi=300)
del axarr, fig

# fig = plt.figure()
# ax = fig.add_subplot(111)
# bins = np.linspace(-1000, 1000, 50)
# for s in slices:
#     mask_1 = cor_czen > slices[s][0]
#     mask_2 = cor_czen <= slices[s][1]
#     czen_mask = np.logical_and(mask_1, mask_2)
#     ax.hist(
#         cor_ca_z[czen_mask], bins=bins, histtype='step',
#         weights=which_weights[czen_mask],
#         label=f"({slices[s][0]}, {slices[s][1]}]"
#     )
# ax.legend()
# ax.set_yscale('log')
# ax.set_xlabel('Closest Approach Depth')
# ax.set_ylabel('Weighted Events')
# fig.tight_layout()
# fig.savefig(f'./plots/ca_z_hist.png')
# del ax, fig


# # ######
# # #### Energy Estimator vs True E
# # ######

# fig2 = plt.figure()
# ax2 = fig2.add_subplot(111)
# bins = [np.linspace(4,10,60), e_plot_bins]
# vals, x_edges, y_edges, im = ax2.hist2d(
#     np.log10(cor_energy), np.log10(cor_energy_est), 
#     weights=which_weights,
#     bins=bins, cmap=my_map, norm=my_norm
# )

# sim_cbar = plt.colorbar(im, ax=ax2, label='unweighted events')
# use_weights = True
# im.set_clim(lims)
# sim_cbar.ax.set_ylabel('weighted events', rotation=270)
    
# ax2.set_xlabel('Primary Energy log10(GeV)')
# ax2.set_ylabel(cfg_file['variables'][which_var]['label'])
# fig2.tight_layout()
# fig2.savefig(f'./plots/{which_var}_vs_energy.png')
# del ax2, fig2, im

######
#### Energy Estimator vs cos(zen)
######

def make_depth_mask(which_slice, depths):
    # for making depth slices
    if which_slice is 'shallow':
        hi = 0
        lo = 2000
    elif which_slice is 'dust':
        hi = -150
        lo = 0
    elif which_slice is 'deep':
        hi = -2000
        lo = -150
    elif which_slice is 'all':
        hi = -100000
        lo = 100000
    hi_mask = depths >= hi
    lo_mask = depths < lo
    return np.logical_and(hi_mask, lo_mask)

depth_slices = {
    0: "all",
    1: "deep",
    2: "dust",
    3: "shallow"
}

rlogl_mask = cor_rlogl < 8.5
speed_mask = cor_speed > 0.28
# speed_mask = np.ones_like(cor_speed)

fig, (ax, ax2) = plt.subplots(1,2, figsize=(7,5))
ax.hist(
    cor_rlogl,
    bins=np.linspace(0,20),
    weights=which_weights
)
ax.set_yscale('log')
ax.set_xlabel("SplineMPE RLOGL")
ax.set_ylabel("Events/Yr")

ax2.hist(
    cor_speed,
    bins=np.linspace(0,0.5, 50),
    weights=which_weights
)
ax2.set_yscale('log')
ax2.set_xlabel('LineFit Speed')
ax2.set_ylabel('Events/Yr')
fig.tight_layout()
fig.savefig(f"./plots/splinempe_rlogl.png", dpi=300)
del fig, ax, ax2

xvals_zenith = np.radians(np.linspace(0,180,181))
func = np.vectorize(helper.get_lognpecut_by_zenith_9yr)
yvals_lognpe = func(xvals_zenith, which_var)

print(f"Summed Weights {np.sum(which_weights)}")
print(np.any(np.isnan(cor_energy_est)))

energyest_vs_czen_bins = [ np.linspace(-1, 1, 20), e_plot_bins ]
fig, axarr =  plt.subplots(2,2, figsize=(8,8))
axarr = axarr.flatten()
for i, d in enumerate(depth_slices.keys()):
    this_depth_mask = make_depth_mask(depth_slices[d], cor_ca_z)
    total_mask = np.logical_and(this_depth_mask, speed_mask)
    total_mask = this_depth_mask
    a, b, c, im = axarr[i].hist2d(
        cor_czen[total_mask], np.log10(cor_energy_est)[total_mask],
        weights=which_weights[total_mask],
        bins=energyest_vs_czen_bins, cmap=my_map, norm=my_norm, cmin=1
    )
    axarr[i].plot(np.cos(xvals_zenith), yvals_lognpe, '-', color='black')
    cbar = plt.colorbar(im, ax=axarr[i], label='Evts / Year')
    im.set_clim(lims)
    axarr[i].set_xlabel(r"SplineMPE cos$\theta$")
    axarr[i].set_ylabel(cfg_file['variables'][which_var]['label'])
    axarr[i].set_title(f"COG in {depth_slices[d]} ice")
fig.tight_layout()
fig.savefig(f"./plots/{which_var}_vs_czen.png", dpi=300)
del fig, axarr


# fig, ax = plt.subplots(1,1, figsize=(7,5))
# a, b, c, im1 = ax.hist2d(
#     cor_czen, np.log10(cor_energy_est), 
#     weights=which_weights,
#     bins = energyest_vs_czen_bins, cmap=my_map, norm=my_norm
# )
# cbar = plt.colorbar(im1, ax=ax, label='Evts / Year')
# im1.set_clim(lims)
# ax.set_xlabel(r"SplineMPE cos$\theta$")
# ax.set_ylabel(cfg_file['variables'][which_var]['label'])
# fig.tight_layout()
# fig.savefig(f"./plots/{which_var}_vs_czen.png", dpi=300)
# del fig, ax