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

readout_file = "/disk20/users/brian/IceCube/juliet/mu_high_energy_merged_1000files.hdf5"

the_f = tables.open_file(readout_file)
charge = the_f.get_node("/Homogenized_QTot").col('value')
errata = the_f.get_node("/LenCalErrata").col('value')
windows = the_f.get_node("/LenSatWindows").col('value')
energy = the_f.get_node("/I3JulietPrimaryParticle").col('energy')
qmin = charge > 2E4

charge_bins = [ np.logspace(3, 8, 16), np.linspace(0, 500) ]
energy_bins = [ np.logspace(4, 12, 16), np.linspace(0, 500) ]
my_map = plt.cm.plasma
norm = colors.LogNorm()

# calibration errata
fig, ((ax1, ax2), (ax3,ax4)) = plt.subplots(2,2, figsize=(12,10))
a, b, c, im = ax1.hist2d(
    charge[qmin], errata[qmin],
    bins = charge_bins, cmap=my_map, norm=norm, cmin=1
)
cbar = plt.colorbar(im, ax=ax1, label='Unweighted Counts')
ax1.set_xlabel("Charge / PE")
ax1.set_ylabel("Num DOMs with Errata")
ax1.set_xscale('log')

a, b, c, im2 = ax2.hist2d(
    energy[qmin], errata[qmin],
    bins = energy_bins, cmap=my_map, norm=norm, cmin=1
)
cbar2 = plt.colorbar(im2, ax=ax2, label='Unweighted Counts')
ax2.set_xlabel("Energy at Detector / GeV")
ax2.set_ylabel("Num DOMs with Errata")
ax2.set_xscale('log')

# saturation windows
a, b, c, im3 = ax3.hist2d(
    charge[qmin], windows[qmin],
    bins = charge_bins, cmap=my_map, norm=norm, cmin=1
)
cbar = plt.colorbar(im, ax=ax1, label='Unweighted Counts')
ax1.set_xlabel("Charge / PE")
ax1.set_ylabel("Num DOMs with Saturation Windows")
ax1.set_xscale('log')

a, b, c, im2 = ax2.hist2d(
    energy[qmin], windows[qmin],
    bins = energy_bins, cmap=my_map, norm=norm, cmin=1
)
cbar2 = plt.colorbar(im2, ax=ax2, label='Unweighted Counts')
ax2.set_xlabel("Energy at Detector / GeV")
ax2.set_ylabel("Num DOMs with Saturation Windows")
ax2.set_xscale('log')



fig.tight_layout()
fig.savefig("cal_issues_vs_charge.png")


the_f.close()