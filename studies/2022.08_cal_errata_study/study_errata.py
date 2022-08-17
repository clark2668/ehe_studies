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
lims = [1E-8, 1E-2]


cfg_file = '../2022.07_cuts/config.yaml'
cfg_file = yaml.safe_load(open(cfg_file))

readout_file = cfg_file['juliet']['mu']['high_energy']['file_cteq5']

livetime = 60*60*24*365

the_f = tables.open_file(readout_file)
charge = the_f.get_node("/Homogenized_QTot").col('value')
errata = the_f.get_node("/LenCalErrata").col('value')
windows = the_f.get_node("/LenSatWindows").col('value')
energy = the_f.get_node("/I3JulietPrimaryParticle").col('energy')
weight_dict, prop_matrix, evts_per_file = weighting.get_juliet_weightdict_and_propmatrix(the_f)
n_gen = cfg_file['juliet']['mu']['high_energy']['n_files'] * evts_per_file
weights = weighting.calc_juliet_weight_from_weight_dict_and_prop_matrix(
    weight_dict=weight_dict, prop_matrix=prop_matrix, 
    flux=gzk_partial, n_gen=n_gen, livetime=livetime
)

qmin = charge > 2E4

charge_bins = [ np.logspace(3, 8, 16), np.linspace(0, 500) ]
energy_bins = [ np.logspace(4, 12, 16), np.linspace(0, 500) ]
my_map = plt.cm.plasma
norm = colors.LogNorm()

# calibration errata
fig, ((ax1, ax2), (ax3,ax4)) = plt.subplots(2,2, figsize=(12,10))
a, b, c, im = ax1.hist2d(
    charge[qmin], errata[qmin], weights=weights[qmin],
    bins = charge_bins, cmap=my_map, norm=norm
)
cbar = plt.colorbar(im, ax=ax1, label='Evts / Year')
im.set_clim(lims)
ax1.set_xlabel("Charge / PE")
ax1.set_ylabel("Num DOMs with Errata")
ax1.set_xscale('log')

a, b, c, im2 = ax2.hist2d(
    energy[qmin], errata[qmin], weights=weights[qmin],
    bins = energy_bins, cmap=my_map, norm=norm
)
cbar2 = plt.colorbar(im2, ax=ax2, label='Evts / Year')
im2.set_clim(lims)
ax2.set_xlabel("Energy at Detector / GeV")
ax2.set_ylabel("Num DOMs with Errata")
ax2.set_xscale('log')

# saturation windows
a, b, c, im3 = ax3.hist2d(
    charge[qmin], windows[qmin], weights=weights[qmin],
    bins = charge_bins, cmap=my_map, norm=norm
)
cbar = plt.colorbar(im3, ax=ax3, label='Evts / Year')
im3.set_clim(lims)
ax3.set_xlabel("Charge / PE")
ax3.set_ylabel("Num DOMs with Saturation Windows")
ax3.set_xscale('log')

a, b, c, im4 = ax4.hist2d(
    energy[qmin], windows[qmin], weights=weights[qmin],
    bins = energy_bins, cmap=my_map, norm=norm
)
cbar2 = plt.colorbar(im4, ax=ax4, label='Evts / Year')
im4.set_clim(lims)
ax4.set_xlabel("Energy at Detector / GeV")
ax4.set_ylabel("Num DOMs with Saturation Windows")
ax4.set_xscale('log')

fig.suptitle(r"$\mu$, Ahlers GZK", size=20)

fig.tight_layout()
fig.savefig("cal_issues_vs_charge.png", dpi=300)


the_f.close()