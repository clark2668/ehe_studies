import h5py
import astropy
from astropy.time import Time
from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np
import argparse

hqtot_key = 'HomogenizedQTot'
hqtot_key = 'HomogenziedQtot_SplitInIcePulses'
# hqtot_key = 'HomogenziedQtot_SplitInIcePulses_ExcludeCalibrationErrata'

#https://wiki.icecube.wisc.edu/index.php/Standard_Candle#Frequently_asked_questions
filter_settings = [1, 2, 3, 4, 5, 6]
filter_brightness = [1.4, 2.92, 8.4, 26, 37.2, 100] 

for f in filter_settings:
	f = h5py.File()

# start the list
f = h5py.File(args.input_files[0], 'r')
charges = f[hqtot_key]['value']
times = f['I3EventHeader']['time_start_mjd']
events = f['I3EventHeader']['Event']
portia_charges = f['EHEPortiaEventSummarySRT']['bestNPEbtw']
print(f['I3EventHeader'].dtype.names)
f.close()

# loop over the remaining files
for temp_f in args.input_files[1:]:
	print("File {}".format(temp_f))
	f = h5py.File(temp_f,'r')
	temp_charges = f[hqtot_key]['value']
	temp_times = f['I3EventHeader']['time_start_mjd']
	temp_portia_charges = f['EHEPortiaEventSummarySRT']['bestNPEbtw']
	temp_events = f['I3EventHeader']['Event']
	charges = np.concatenate((charges, temp_charges))
	times = np.concatenate((times,temp_times))
	portia_charges = np.concatenate((portia_charges,temp_portia_charges))
	events = np.concatenate((events,temp_events))

# mask = portia_charges > 4E4
# mask = portia_charges > 0
# portia_charges = portia_charges[mask]
# charges = charges[mask]


## quick study of the "flatline" events (evts with portia Q by not Homog Q)
mask = charges < 10000
portia_charges = portia_charges[mask]
charges = charges[mask]
events = events[mask]
for p, h, e in zip(portia_charges, charges, events):
	if p > 20000:
		print("Ev {}, Portia {:.2f}, Hqtot {:.2f}, H/P {:.2f}".format(e, p, h, h/p))

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


