import numpy as np
import h5py as h5py
import matplotlib.pyplot as plt
import sys

print(sys.argv)
use_fadc=bool(int(sys.argv[1]))
use_atwd=bool(int(sys.argv[2]))
print('use fadc {}'.format(use_fadc))
print('use atwd {}'.format(use_atwd))
hese_pulses = 'SplitInIcePulses'

print("HESE Pulses: {}, FADC: {}".format(hese_pulses, use_fadc))
print("----------------------------")

file_in = h5py.File('comparison_overlap_{}_FADC_{}_ATWD_{}.hdf5'.format(hese_pulses, use_fadc, use_atwd), 'r')
data_overlap = file_in['data_overlap']
hese_npe = np.asarray(data_overlap['hese_overlap'])
ehe_npe = np.asarray(data_overlap['ehe_overlap'])
strings = np.asarray(data_overlap['strings'])
doms = np.asarray(data_overlap['doms'])

baselines_in = h5py.File('compare_baselines.hdf5', 'r')
data_baselines = baselines_in['baselines']
strings_baselines = data_baselines['strings']
doms_baselines = data_baselines['doms']
portia_fadc_baseline = np.asarray(data_baselines['portia_fadc'])
portia_atwd_baseline = np.asarray(data_baselines['portia_atwd'])
baselines_portia_fadc = {}
baselines_portia_atwd = {}
for i, base in enumerate(portia_fadc_baseline):
	baselines_portia_fadc[(int(strings_baselines[i]), int(doms_baselines[i]))] = portia_fadc_baseline[i]
for i, base in enumerate(portia_atwd_baseline):
	baselines_portia_atwd[(int(strings_baselines[i]), int(doms_baselines[i]))] = portia_atwd_baseline[i]

# calculate differences
diff = []
big_ehe = []
big_hese = []
for this_ehe, this_hese, string, dom in zip(ehe_npe, hese_npe, strings, doms):
	rel_diff = (this_hese - this_ehe)/this_ehe
	# print("({}, {}): EHE {:.2f}, HESE {:.2f}, rel dif {:.2f}".format(string, dom, this_ehe, this_hese, rel_diff))	
	if rel_diff < -0.9:
		# print("({}, {}): EHE {:.2f}, HESE {:.2f}, rel dif {:.2f}".format(string, dom, this_ehe, this_hese, rel_diff))
		big_ehe.append(this_ehe)
		big_hese.append(this_hese)
	diff.append(rel_diff)
diff = np.asarray(diff)
big_ehe = np.asarray(big_ehe)
big_hese = np.asarray(big_hese)

plot_baseline_stuff = True
if plot_baseline_stuff:
	fig, axs = plt.subplots(1,1,figsize=(5,5))
	axs.hist(portia_fadc_baseline/1e-12, bins=50, alpha=0.5, label='N={}'.format(len(diff)))#, histtype='step')
	# axs.set_yscale('log')
	axs.set_xlabel('Baselines')
	axs.set_ylabel('Number of FADCs')
	# axs.set_xlim([-1.0, 1.0])
	# axs.set_title('FADC = {}, ATWD = {}'.format(use_fadc, use_atwd))
	# axs.set_ylim([0,60])
	# axs.plot([0,250],[0,250],'--', label='1-1 line')
	axs.legend()
	plt.tight_layout()
	fig.savefig('fadc_baselines.png', dpi=300)
	del fig, axs

	these_fadcs = []
	these_atwds = []
	for string, dom in zip(strings, doms):
		this_fadc = baselines_portia_fadc[(int(string), int(dom))]
		this_atwd = baselines_portia_atwd[(int(string), int(dom))]
		these_fadcs.append(this_fadc)
		these_atwds.append(this_atwd)
		print("om {}, string {}, fadc {}, atwd {}, fadc/atwd {}".format(string, dom, this_fadc, this_atwd, this_fadc/this_atwd ))
	these_fadcs = np.asarray(these_fadcs)

	fig, axs = plt.subplots(1,1,figsize=(5,5))
	# axs.plot(these_fadcs, diff, 'o', alpha=0.5, label='N={}'.format(len(diff)))#, histtype='step')
	# axs.plot(these_fadcs, ehe_npe/hese_npe, 'o', alpha=0.5, label='N={}'.format(len(diff)))
	# axs.plot(portia_atwd_baseline/1E-12, portia_fadc_baseline/1E-12, 'o', alpha=0.5, label='N={}'.format(len(portia_atwd_baseline)))
	axs.plot(these_fadcs/these_atwds, ehe_npe/hese_npe, 'o', alpha=0.5, label='N={}'.format(len(these_fadcs)))
	axs.set_xlabel('FADC/ATWD Baseline')
	axs.set_ylabel('EHE/HESE NPE')


	# axs.set_yscale('log')
	# axs.set_xlabel('Portia FADC Baseline')
	# axs.set_ylabel('Relative Difference')
	# axs.set_ylabel('EHE/HESE')
	# axs.set_xlim([-0.5, 0.5])
	# axs.set_ylim([-1.5, 1.5])
	# axs.set_title('FADC = {}, ATWD = {}'.format(use_fadc, use_atwd))
	# axs.set_ylim([0,60])
	# axs.plot([0,250],[0,250],'--', label='1-1 line')
	axs.legend()
	plt.tight_layout()
	fig.savefig('fadc_to_atwd__vs__ehe_to_hese.png', dpi=300)
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
	axs.set_title('FADC = {}, ATWD = {}'.format(use_fadc, use_atwd))
	axs.legend()
	axs.plot(big_ehe, big_hese, 'o', alpha=0.5, color='red')

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
	fig.savefig('ehenpe_vs_hesenpe_{}_FADC_{}_ATWD_{}_log.png'.format(hese_pulses, use_fadc, use_atwd), dpi=300)
	plt.close(fig)
	del fig, axs


	fig, axs = plt.subplots(1,1,figsize=(5,5))
	axs.hist(diff, bins=50, alpha=0.5, label='N={}'.format(len(diff)))#, histtype='step')
	# axs.set_yscale('log')
	axs.set_xlabel('(HESE - EHE)/EHE  [HESE={}]'.format(hese_pulses))
	axs.set_ylabel('Number of DOMs')
	axs.set_xlim([-1.0, 1.0])
	axs.set_title('FADC = {}, ATWD = {}'.format(use_fadc, use_atwd))
	axs.set_ylim([0,60])
	# axs.plot([0,250],[0,250],'--', label='1-1 line')
	axs.legend()
	plt.tight_layout()
	fig.savefig('ehe_vs_hese_rel_diff_{}_FADC_{}_ATWD_{}.png'.format(hese_pulses, use_fadc, use_atwd), dpi=300)

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
	axs.set_title('FADC = {}, ATWD = {}'.format(use_fadc, use_atwd))
	axs.set_xlabel('NPE')
	axs.set_ylabel('Number of DOMs')
	# axs.set_yscale('log')
	axs.set_xscale('log')
	axs.legend()
	plt.tight_layout()
	fig.savefig('ehe_vs_hese_missing_{}_FADC_{}_ATWD_{}.png'.format(hese_pulses, use_fadc, use_atwd), dpi=300)

