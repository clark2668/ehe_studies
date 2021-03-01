import numpy as np
import h5py as h5py
import matplotlib.pyplot as plt
import sys

print(sys.argv)
use_fadc=bool(int(sys.argv[1]))
use_atwd=bool(int(sys.argv[2]))
hese_pulses = 'SplitInIcePulses'
beacon_fadc=bool(int(sys.argv[3])) #True
noise_cut=bool(int(sys.argv[4])) #True

print("HESE Pulses: {}, Use FADC {}, Use ATWD {}, Beacon FADC {}, Noise Cut {}".format(hese_pulses, use_fadc, use_atwd, beacon_fadc, noise_cut))
print("----------------------------")

# file_in = h5py.File('comparison_overlap_{}_FADC_{}_ATWD_{}_forcezero_{}.hdf5'.format(hese_pulses, use_fadc, use_atwd, force_zero), 'r')
# file_in = h5py.File('comparison_overlap_{}_FADC_{}_ATWD_{}_use_beacon_fadc.hdf5'.format(hese_pulses, use_fadc, use_atwd), 'r')
file_in = h5py.File('comparison_overlap_HESE_{}_FADC_{}_ATWD_{}_FADCBeacon_{}_NoiseCut_{}.hdf5'.format(hese_pulses, use_fadc, use_atwd, beacon_fadc, noise_cut), 'r')
data_overlap = file_in['data_overlap']
hese_npe = np.asarray(data_overlap['hese_overlap'])
ehe_npe = np.asarray(data_overlap['ehe_overlap'])
strings = np.asarray(data_overlap['strings'])
doms = np.asarray(data_overlap['doms'])

# calculate differences
diff = []
big_ehe = []
big_hese = []
ratio = []
for this_ehe, this_hese, string, dom in zip(ehe_npe, hese_npe, strings, doms):
	rel_diff = (this_hese - this_ehe)/this_ehe
	# print("({}, {}): EHE {:.2f}, HESE {:.2f}, rel dif {:.2f}".format(string, dom, this_ehe, this_hese, rel_diff))	
	if rel_diff < -0.7:
		print("({}, {}): EHE {:.2f}, HESE {:.2f}, rel dif {:.2f}".format(string, dom, this_ehe, this_hese, rel_diff))
		big_ehe.append(this_ehe)
		big_hese.append(this_hese)
	diff.append(rel_diff)
	ratio.append(this_ehe/this_hese)
diff = np.asarray(diff)
big_ehe = np.asarray(big_ehe)
big_hese = np.asarray(big_hese)
ratio = np.asarray(ratio)

plot_baseline_stuff = True
if plot_baseline_stuff:

	baselines_in = h5py.File('compare_baselines.hdf5', 'r')
	data_baselines = baselines_in['baselines']
	strings_baselines = data_baselines['strings']
	doms_baselines = data_baselines['doms']

	portia_fadc_baseline = np.asarray(data_baselines['portia_fadc'])
	beacon_fadc_baseline = np.asarray(data_baselines['beacon_fadc'])

	fig, axs = plt.subplots(1,1,figsize=(5,5))
	axs.plot(portia_fadc_baseline, beacon_fadc_baseline, 'o', alpha=0.5)
	axs.set_xlabel('Portia FADC Baseline [mV]')
	axs.set_ylabel('Beacon FADC Base [mV]')
	axs.set_xlim([-0.15, 0.15])
	axs.set_ylim([-0.15, 0.15])
	axs.set_aspect('equal')
	plt.tight_layout()
	fig.savefig('fadc_baseline_comparison.png', dpi=300)
	del fig, axs


	# and a histogram
	fig, axs = plt.subplots(1,1,figsize=(5,5))
	bins = np.linspace(-0.15,0.15,100)
	axs.hist(portia_fadc_baseline, bins=bins, alpha=0.5, label='Portia')#, histtype='step')
	axs.hist(beacon_fadc_baseline, bins=bins, alpha=0.5, label='Beacon')#, histtype='step')
	axs.set_xlabel('FADC Baseline [mV]')
	axs.set_ylabel('Number of DOMs')
	axs.set_yscale('log')
	axs.set_xlim([-0.15, 0.15])
	# axs.set_title('FADC = {}, ATWD = {}'.format(use_fadc, use_atwd))
	# axs.plot([0,250],[0,250],'--', label='1-1 line')
	axs.legend()
	plt.tight_layout()
	fig.savefig('fadc_baseline_comparison_hist.png', dpi=300)
	del fig, axs


plot_ehe_and_hese_comparison=False
if plot_ehe_and_hese_comparison:

	# make a plot of EHE Q vs HESE Q
	fig, axs = plt.subplots(1,1,figsize=(5,5))
	print("Number of mutual DOMs is {}".format(len(ehe_npe)))
	axs.plot(ehe_npe, hese_npe,'o', alpha=0.5, label='N={}'.format(len(ehe_npe)))
	axs.set_xlabel('EHE NPE')
	axs.set_ylabel('HESE NPE ({})'.format(hese_pulses))
	axs.plot([0,300],[0,300],'--', label='1-1 line')
	# axs.set_title('FADC = {}, ATWD = {}, FADC Beacon Baselines {}, Noise Cut'.format(use_fadc, use_atwd, beacon_fadc, noise_cut))
	axs.set_title('Use FADC = {}, Use ATWD = {}, \nUse FADC Beacon Baselines={}, \n Use FADC Noise Cut={}'.format(use_fadc, use_atwd, beacon_fadc, noise_cut), size=12)
	axs.legend()
	# axs.plot(big_ehe, big_hese, 'o', alpha=0.5, color='red')

	axs.set_xlim([-10, 300])
	axs.set_ylim([-10, 300])
	axs.legend()
	# axs.set_aspect('equal')
	# plt.tight_layout()
	# fig.savefig('ehenpe_vs_hesenpe_{}_FADC_{}_ATWD_{}.png'.format(hese_pulses, use_fadc, use_atwd), dpi=300)

	axs.set_yscale('log')
	axs.set_xscale('log')
	axs.set_xlim([0.1, 300])
	axs.set_ylim([0.1, 300])
	axs.set_aspect('equal')
	plt.tight_layout()
	# fig.savefig('ehenpe_vs_hesenpe_{}_FADC_{}_ATWD_{}__forcezero_{}log.png'.format(hese_pulses, use_fadc, use_atwd, force_zero), dpi=300)
	# fig.savefig('ehenpe_vs_hesenpe_{}_FADC_{}_ATWD_{}_use_beacon_fadc_log.png'.format(hese_pulses, use_fadc, use_atwd), dpi=300)
	fig.savefig('ehenpe_vs_hesenpe_HESE_{}_FADC_{}_ATWD_{}_FADCBeacon_{}_NoiseCut_{}.png'.format(hese_pulses, use_fadc, use_atwd, beacon_fadc, noise_cut), dpi=300)
	plt.close(fig)
	del fig, axs

	bins = np.linspace(0,25,26)
	fig, axs = plt.subplots(1,1,figsize=(5,5))
	axs.hist(ratio, bins=bins, alpha=0.5, label='N={}'.format(len(ratio)))#, histtype='step')
	# axs.set_yscale('log')
	axs.set_xlabel('EHE/HESE  [HESE={}]'.format(hese_pulses))
	axs.set_ylabel('Number of DOMs')
	# axs.set_xlim([-1.0, 1.0])
	axs.set_title('Use FADC = {}, Use ATWD = {}, \nUse FADC Beacon Baselines={}, \n Use FADC Noise Cut={}'.format(use_fadc, use_atwd, beacon_fadc, noise_cut), size=12)
	# axs.set_ylim([0,60])
	# axs.plot([0,250],[0,250],'--', label='1-1 line')
	axs.set_yscale('log')
	axs.legend()
	plt.tight_layout()
	# fig.savefig('ehe_over_hese_{}_FADC_{}_ATWD_{}_forcezero_{}.png'.format(hese_pulses, use_fadc, use_atwd, force_zero), dpi=300)
	fig.savefig('ehenpe_over_hese_HESE_{}_FADC_{}_ATWD_{}_FADCBeacon_{}_NoiseCut_{}.png'.format(hese_pulses, use_fadc, use_atwd, beacon_fadc, noise_cut), dpi=300)
	del fig, axs

	##############################
	##############################

	fig, axs = plt.subplots(1,1,figsize=(5,5))
	axs.hist(diff, bins=50, alpha=0.5, label='N={}'.format(len(diff)))#, histtype='step')
	# axs.set_yscale('log')
	axs.set_xlabel('(HESE - EHE)/EHE  [HESE={}]'.format(hese_pulses))
	axs.set_ylabel('Number of DOMs')
	axs.set_xlim([-1.0, 1.0])
	axs.set_title('Use FADC = {}, Use ATWD = {}, \nUse FADC Beacon Baselines={}, \n Use FADC Noise Cut={}'.format(use_fadc, use_atwd, beacon_fadc, noise_cut), size=12)
	axs.set_ylim([0,60])
	# axs.plot([0,250],[0,250],'--', label='1-1 line')
	axs.legend()
	plt.tight_layout()
	# fig.savefig('ehe_vs_hese_rel_diff_{}_FADC_{}_ATWD_{}_forcezero_{}.png'.format(hese_pulses, use_fadc, use_atwd, force_zero), dpi=300)
	fig.savefig('ehe_vs_hese_reldiff_HESE_{}_FADC_{}_ATWD_{}_FADCBeacon_{}_NoiseCut_{}.png'.format(hese_pulses, use_fadc, use_atwd, beacon_fadc, noise_cut), dpi=300)

	##############################
	##############################	

	data_exclusion = file_in['data_exclusion']
	missing_in_hese = np.asarray(data_exclusion['missing_in_hese'])
	missing_in_ehe = np.asarray(data_exclusion['missing_in_ehe'])

	fig, axs = plt.subplots(1,1,figsize=(5,5))
	# bins = range(-10,1000+10,10)
	bins = np.logspace(np.log10(0.1), np.log10(1e3), 50)
	axs.hist(missing_in_hese, bins=bins, alpha=0.5, label='Missing in HESE')#, histtype='step')
	axs.hist(missing_in_ehe, bins=bins, alpha=0.5, label='Missing in EHE ({})'.format(hese_pulses))#, histtype='step')
	axs.set_title('Use FADC = {}, Use ATWD = {}, \nUse FADC Beacon Baselines={}, \n Use FADC Noise Cut={}'.format(use_fadc, use_atwd, beacon_fadc, noise_cut), size=12)
	axs.set_xlabel('NPE')
	axs.set_ylabel('Number of DOMs')
	# axs.set_yscale('log')
	axs.set_xscale('log')
	axs.legend()
	plt.tight_layout()
	# fig.savefig('ehe_vs_hese_missing_{}_FADC_{}_ATWD_{}_forcezero_{}.png'.format(hese_pulses, use_fadc, use_atwd, force_zero), dpi=300)
	fig.savefig('ehe_vs_hese_missing_reldiff_HESE_{}_FADC_{}_ATWD_{}_FADCBeacon_{}_NoiseCut_{}.png'.format(hese_pulses, use_fadc, use_atwd, beacon_fadc, noise_cut), dpi=300)

