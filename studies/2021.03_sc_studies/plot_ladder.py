import h5py
import astropy
from astropy.time import Time
from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np
import argparse

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-i", type=str,
	dest="file_location", default='/Users/brianclark/Documents/work/IceCube/ehe/step4_files',
	help="path to the hdf5 file locations with HQtot etc saved out",
	)
parser.add_argument("-y", type=int,
	dest="year", default=2015,
	help="Which year of standard candle, e.g. 2015"
	)
parser.add_argument("-c", type=int,
	dest="candle", default=2,
	help="Which candle, e.g. 2"
	)
parser.add_argument("-o", type=str,
	dest="option", required=True,
	help="which option, standard or magsix"
	)
parser.add_argument("-p", type=str,
	dest="pulses", required=True,
	help="which pulses to choose, hqtot or portia"
	)
args = parser.parse_args()

option = args.option
if option=='standard':
	hqtot_key = 'HomogenziedQtot_SplitInIcePulses'
	portia_key = 'EHEPortiaEventSummarySRT'
elif option=='magsix':
	hqtot_key = 'HomogenizedQTot_DeepMagSix'
	portia_key = 'PortiaEventSummarySRT_DeepMagSix'
pulses = args.pulses

#https://wiki.icecube.wisc.edu/index.php/Standard_Candle#Frequently_asked_questions
filter_settings = [1, 2, 3, 4, 5, 6]
filter_brightness = [1.4, 2.92, 8.4, 26, 37.2, 100]
# filter_brightness = [1.4, 3.28, 8.90, 26.8, 37.7, 100]
brightness_dict = {}
for i, j in zip(filter_settings, filter_brightness):
	brightness_dict[i] = j * 1E-2

def get_predicted_brightness(original_brightness, filter_setting):
	return original_brightness * brightness_dict[filter_setting]

setting_charge_dict = {}

for f in filter_settings:
	file = h5py.File(args.file_location + '/' + f'y{args.year}_c{args.candle}_f{f}_hqtot_splitinnicepulses.hdf5', 'r')
	if pulses=='hqtot':
		charges = file[hqtot_key]['value']
	elif pulses=='portia':
		charges = file[portia_key]['bestNPEbtw']
	file.close()
	new_charges = []
	for charge in charges:
		if charge > 1E3 and f<4:
			new_charges.append(charge)
		elif f>3 and charge > 1E4: # manually remove "flatline" events (this check can probably go away later)
			new_charges.append(charge)
	setting_charge_dict[f] = new_charges

# start with filter setting 1, which is the "bottom" rung of the brightness ladder
bins = np.linspace(4.75,4.95,50)
if option == 'magsix':
	bins = np.linspace(3.5, 4.15, 50)

fig, axs = plt.subplots(1,1,figsize=(5,5))
filter_1_brightness = np.average(setting_charge_dict[1])
original_brightness = filter_1_brightness/brightness_dict[1]
print('Original brightness {}'.format(original_brightness))
axs.hist(np.log10(setting_charge_dict[1]), bins=bins, alpha=0.5, label='Filter Setting 1')
axs.vlines(np.log10(filter_1_brightness), 0.1, 1E3, label='Mean')
axs.set_yscale('log')
axs.set_ylabel('Number of Events')
axs.set_xlabel('Charge')
axs.legend()
axs.set_title('{}, {}'.format(option, pulses))
plt.tight_layout()
fig.savefig('ladder_plots/hist_of_q_brightness_1_{}_{}.png'.format(option, pulses), dpi=300)
del fig, axs

# calculate predicted brightness assuming no saturatione effects
predicted_brightness_dict = {}
for i in filter_settings:
	predicted_brightness_dict[i] = get_predicted_brightness(original_brightness, i)

colors = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5']


# now, plot the first 3
bins = np.linspace(3.5,6.5,150)

fig, axs = plt.subplots(1,1,figsize=(5,5))
for i, f in enumerate(filter_settings[:6]):
	obsv_mean = np.average(setting_charge_dict[f])
	pred_mean = predicted_brightness_dict[f]
	print("Setting {}, Pred {:.2f}, Obsv {:.2f}, (O-P)/P {:.3f}".format(f, 
		pred_mean, obsv_mean, (obsv_mean-pred_mean)/pred_mean))
	a = axs.hist(np.log10(setting_charge_dict[f]), bins=bins, color=colors[i], alpha=0.5, label='Filter {}'.format(f))
	axs.vlines(np.log10(predicted_brightness_dict[f]), 0.1, 1E3, color=colors[i])
axs.set_yscale('log')
axs.set_ylabel('Number of Events')
axs.set_xlabel('Charge')
axs.set_title('{}, {}'.format(option, pulses))
axs.legend()
plt.tight_layout()
fig.savefig('ladder_plots/hist_of_q_brightness_multi_{}_{}.png'.format(option, pulses), dpi=300)
del fig, axs



