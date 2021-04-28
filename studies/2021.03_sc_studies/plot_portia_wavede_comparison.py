import h5py
import astropy
from astropy.time import Time
from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.patches import Rectangle
from matplotlib.dates import (YEARLY, DateFormatter, rrulewrapper, RRuleLocator, drange)
import numpy as np
import argparse

import tools

parser = argparse.ArgumentParser()
parser.add_argument("-i", type=str, nargs='+',
	dest="input_files",
	help="full path to the input file",
	required=True)
args = parser.parse_args()

# start the list
f = h5py.File(args.input_files[0], 'r')
charges = f['HomogenizedQTot']['value']
times = f['I3EventHeader']['time_start_mjd']
portia_charges = f['EHEPortiaEventSummarySRT']['bestNPEbtw']
f.close()

# loop over the remaining files
for temp_f in args.input_files[1:]:
	print("File {}".format(temp_f))
	f = h5py.File(temp_f,'r')
	temp_charges = f['HomogenizedQTot']['value']
	temp_times = f['I3EventHeader']['time_start_mjd']
	temp_portia_charges = f['EHEPortiaEventSummarySRT']['bestNPEbtw']
	charges = np.concatenate((charges, temp_charges))
	times = np.concatenate((times,temp_times))
	portia_charges = np.concatenate((portia_charges,temp_portia_charges))


# plot the charge in histogram
# NB: by construction, this probably excludes balloon DOMs, which might be interesting...
# fig, axs = plt.subplots(1,1,figsize=(5,5))
# bins = np.linspace(0.5E5, 1.2E5, 50)
# axs.hist(charges, bins=bins,alpha=0.5, label='HomogenizedQTot')
# axs.hist(portia_charges, bins=bins, alpha=0.5, label='BestNPEBTW')
# axs.set_yscale('log')
# axs.set_ylabel('Number of Events')
# axs.set_xlabel('Charge')
# plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
# axs.legend()
# plt.tight_layout()
# fig.savefig('hist_of_charge_portia_homogqtot.png', dpi=300)
# del fig, axs

# plt.ticklabel_format(axis="x", style="plain") # reset

# and as a scatter plot

fig, axs = plt.subplots(1,1,figsize=(5,5))
axs.plot(np.log10(portia_charges), np.log10(charges), 'o', alpha=0.5, label='SC Events')
axs.plot([0,1E6],[0,1E6],'--', label='1-1 line')
# axs.set_yscale('log')
axs.set_ylabel('Homogenized Q Tot')
axs.set_xlabel('Portia Best NPE BTW')
# axs.set_xlim([5E4, 80E4])
# axs.set_ylim([5E4, 80E4])
axs.set_aspect('equal')
axs.legend()
# plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
plt.tight_layout()
fig.savefig('portia_vs_homogqtot.png', dpi=300)
del fig, axs
