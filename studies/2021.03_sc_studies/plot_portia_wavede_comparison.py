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

parser = argparse.ArgumentParser()
parser.add_argument("-i", type=str, nargs='+',
	dest="input_files",
	help="full path to the input file",
	required=True)
args = parser.parse_args()

hqtot_key = 'HomogenizedQTot'
hqtot_key = 'HomogenziedQtot_SplitInIcePulses'
hqtot_key = 'HomogenziedQtot_SplitInIcePulses_ExcludeCalibrationErrata'


# start the list
f = h5py.File(args.input_files[0], 'r')
charges = f[hqtot_key]['value']
times = f['I3EventHeader']['time_start_mjd']
portia_charges = f['EHEPortiaEventSummarySRT']['bestNPEbtw']
f.close()

# loop over the remaining files
for temp_f in args.input_files[1:]:
	print("File {}".format(temp_f))
	f = h5py.File(temp_f,'r')
	temp_charges = f[hqtot_key]['value']
	temp_times = f['I3EventHeader']['time_start_mjd']
	temp_portia_charges = f['EHEPortiaEventSummarySRT']['bestNPEbtw']
	charges = np.concatenate((charges, temp_charges))
	times = np.concatenate((times,temp_times))
	portia_charges = np.concatenate((portia_charges,temp_portia_charges))

mask = portia_charges > 4E4
# mask = portia_charges > 0
portia_charges = portia_charges[mask]
charges = charges[mask]

# for p, h in zip(portia_charges[-30:], charges[-30:]):
# 	print("Portia {:.2f}, Hqtot {:.2f}, H/P {:.2f}".format(p, h, h/p))

log_portia_charges = np.log10(portia_charges)
log_charges = np.log10(charges)

do_hist = False
if do_hist:

	# plot the charge in histogram
	# NB: by construction, this probably excludes balloon DOMs, which might be interesting...
	fig, axs = plt.subplots(1,1,figsize=(5,5))
	bins = np.linspace(4.5, 6, 250)
	axs.hist(log_charges, bins=bins,alpha=0.5, label='HomogenizedQTot')
	axs.hist(log_portia_charges, bins=bins, alpha=0.5, label='BestNPEBTW')
	axs.set_yscale('log')
	axs.set_ylabel('Number of Events')
	axs.set_xlabel('Charge')
	axs.set_ylim([0.9,3E3])
	axs.set_title(hqtot_key, fontsize=8)
	plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
	axs.legend()
	plt.tight_layout()
	fig.savefig('hist_of_charge_portia_homogqtot.png', dpi=300)
	del fig, axs

# plt.ticklabel_format(axis="x", style="plain") # reset

do_scatter = True
if do_scatter:

	do_log = True

	# and as a scatter plot

	fig, axs = plt.subplots(1,1,figsize=(5,5))
	if do_log:
		axs.plot(np.log10(portia_charges), np.log10(charges), 'o', alpha=0.5, label='SC Events')
		axs.plot([0,6],[0,6],'--', label='1-1 line')
		axs.set_ylabel(r'log$_{10}$(Homogenized Q Tot) [NPE]')
		axs.set_xlabel(r'log$_{10}$(Portia Best NPE BTW) [NPE]')
		axs.set_xlim([4.5, 6])
		axs.set_ylim([4.5, 6])
	else:
		axs.plot(portia_charges, charges, 'o', alpha=0.5, label='SC Events')
		axs.plot([0,1E6],[0,1E6],'--', label='1-1 line')
		axs.set_ylabel(r'Homogenized Q Tot [NPE]')
		axs.set_xlabel(r'Portia Best NPE BTW [NPE]')
		# axs.set_yscale('log')
		axs.set_xlim([0, 1E6])
		axs.set_ylim([0, 1E6])
	axs.set_title(hqtot_key, fontsize=8)
	axs.set_aspect('equal')
	axs.legend()
	plt.tight_layout()
	if do_log:
		fig.savefig('portia_vs_homogqtot_log.png', dpi=300)
	else:
		fig.savefig('portia_vs_homogqtot.png', dpi=300)
	del fig, axs

