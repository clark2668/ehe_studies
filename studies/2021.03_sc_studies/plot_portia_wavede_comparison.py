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

option = 'standard'
if option=='standard':
	hqtot_key = 'HomogenziedQtot_SplitInIcePulses'
	portia_key = 'EHEPortiaEventSummarySRT'
elif option=='magsix':
	hqtot_key = 'HomogenizedQTot_DeepMagSix'
	portia_key = 'PortiaEventSummarySRT_DeepMagSix'


# start the list
f = h5py.File(args.input_files[0], 'r')
charges = f[hqtot_key]['value']
times = f['I3EventHeader']['time_start_mjd']
events = f['I3EventHeader']['Event']
portia_charges = f[portia_key]['bestNPEbtw']
print(f['I3EventHeader'].dtype.names)
f.close()

# loop over the remaining files
for temp_f in args.input_files[1:]:
	print("File {}".format(temp_f))
	f = h5py.File(temp_f,'r')
	temp_charges = f[hqtot_key]['value']
	temp_times = f['I3EventHeader']['time_start_mjd']
	temp_portia_charges = f[portia_key]['bestNPEbtw']
	temp_events = f['I3EventHeader']['Event']
	charges = np.concatenate((charges, temp_charges))
	times = np.concatenate((times,temp_times))
	portia_charges = np.concatenate((portia_charges,temp_portia_charges))
	events = np.concatenate((events,temp_events))

# mask = portia_charges > 4E4
mask = (portia_charges > 0) & (charges > 0)
portia_charges = portia_charges[mask]
charges = charges[mask]


## quick study of the "flatline" events (evts with portia Q by not Homog Q)
# mask = charges < 10000
# portia_charges = portia_charges[mask]
# charges = charges[mask]
# events = events[mask]
# for p, h, e in zip(portia_charges, charges, events):
# 	if p > 20000:
# 		print("Ev {}, Portia {:.2f}, Hqtot {:.2f}, H/P {:.2f}".format(e, p, h, h/p))



do_hist = True
if do_hist:

	log_portia_charges = np.log10(portia_charges)
	log_charges = np.log10(charges)
	
	# to cleanup the plot (temporary!)
	# put a mask on the "flatline" events where Portia is calculating a bunch of chrage
	# so, put a cut on the true charge (hqtot charge) in the frame to be conservative
	mask = log_charges > 3.5
	log_portia_charges = log_portia_charges[mask]
	log_charges = log_charges[mask]


	# plot the charge in histogram
	fig, axs = plt.subplots(1,1,figsize=(5,5))
	bins = np.linspace(3.5, 6, 300)
	axs.hist(log_charges, bins=bins,alpha=0.5, label='HomogenizedQTot')
	axs.hist(log_portia_charges, bins=bins, alpha=0.5, label='BestNPEBTW')
	axs.set_yscale('log')
	axs.set_ylabel('Number of Events')
	axs.set_xlabel('Charge')
	axs.set_ylim([0.9,3E3])
	axs.set_title(option, fontsize=8)
	plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
	axs.legend()
	plt.tight_layout()
	fig.savefig('hist_of_charge_portia_{}_homogqtot_{}.png'.format(portia_key, hqtot_key), dpi=300)
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
		axs.set_ylabel(r'log$_{10}$(Homogenized Q Tot) [NPE]'+'\n'+'{}'.format(hqtot_key))
		axs.set_xlabel(r'log$_{10}$(Portia Best NPE BTW) [NPE]'+'\n'+'{}'.format(portia_key))
		axs.set_xlim([3.5, 6])
		axs.set_ylim([3.5, 6])
	else:
		axs.plot(portia_charges, charges, 'o', alpha=0.5, label='SC Events')
		axs.plot([0,1E6],[0,1E6],'--', label='1-1 line')
		axs.set_ylabel(r'Homogenized Q Tot [NPE]''\n''{}'.format(hqtot_key))
		axs.set_xlabel(r'Portia Best NPE BTW [NPE]''\n''{}'.format(portia_key))
		# axs.set_yscale('log')
		axs.set_xlim([0, 100000])
		axs.set_ylim([0, 100000])
	axs.set_title(option, fontsize=8)
	axs.set_aspect('equal')
	axs.legend()
	plt.tight_layout()
	if do_log:
		fig.savefig('portia_{}_vs_homogqtot_{}_log.png'.format(portia_key, hqtot_key), dpi=300)
	else:
		fig.savefig('portia_{}_vs_homogqtot_{}.png'.format(portia_key, hqtot_key), dpi=300)
	del fig, axs

