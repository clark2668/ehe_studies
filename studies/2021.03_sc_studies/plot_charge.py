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
f.close()

# loop over the remaining files
for temp_f in args.input_files[1:]:
	print("File {}".format(temp_f))
	f = h5py.File(temp_f,'r')
	temp_charges = f['HomogenizedQTot']['value']
	temp_times = f['I3EventHeader']['time_start_mjd']
	charges = np.concatenate((charges, temp_charges))
	times = np.concatenate((times,temp_times))
	# useful for figuring out how to set time window cuts
	# offset = 5.7038E4
	# for t, q in zip(temp_times, temp_charges):
		# if t > offset + 0.195 and t < offset + 0.22 and q > 5E5:
			# print("Something is up , t {}, q {}".format(t+offset, q))

# # plot the charge in histogram
# fig, axs = plt.subplots(1,1,figsize=(5,5))
# axs.hist(charges, bins=100,alpha=0.5)
# axs.set_yscale('log')
# axs.set_ylabel('Number of Events')
# axs.set_xlabel('HomogenizedQTot')
# plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
# plt.tight_layout()
# fig.savefig('hist_of_charge.png', dpi=300)
# del fig, axs

# plt.ticklabel_format(axis="x", style="plain") # reset

# and as a 2D plot

mask = charges > 0.0 # mask "small" charge events
fig, axs = plt.subplots(1,1,figsize=(7,5))
my_map = plt.cm.plasma
counts, xedges, yedges, im = axs.hist2d(times[mask], charges[mask],
	bins=500,
	cmap = my_map,
	norm=colors.LogNorm(),
	cmin=1
	)
cbar = plt.colorbar(im, ax=axs)
cbar.set_label('Number of Events')
# formatter = DateFormatter('%H:%M:%S')
# axs.xaxis.set_major_formatter(formatter)
axs.xaxis.set_tick_params(rotation=45, labelsize=7)
axs.set_ylabel('HomogenizedQTot [NPE]')
axs.set_xlabel('Time [Hour:Min:Sec]')
axs.grid(which='both')
# axs.set_ylim([0,1e5])
plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
plt.tight_layout()
axs.set_ylim(0E5,8E5)


for i in range(1,7,1):
	start, stop, qmin, qmax = tools.get_start_stop(2015, 2, i)
	qmin = 1
	qmax = 8E5
	rect = Rectangle((start, qmin), (stop-start), (qmax-qmin),
		linewidth=1, edgecolor='r', facecolor='None',linestyle='--')
	axs.add_patch(rect)

fig.savefig('charge_vs_time.png'.format(filter), dpi=300)

# might be helpful later: https://github.com/toej93/thesis_work_daily/blob/master/2021_Mar18_CoordSystemSourceSearch.ipynb


