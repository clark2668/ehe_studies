import numpy as np
import h5py as h5py
import matplotlib.pyplot as plt
import sys
import matplotlib.colors as colors
import statistics

hese_pulses = 'SplitInIcePulses'
print(sys.argv)
use_fadc=bool(int(sys.argv[1]))
use_atwd=bool(int(sys.argv[2]))
beacon_fadc=bool(int(sys.argv[3])) #True
noise_cut=bool(int(sys.argv[4])) #True
causal_qtot=bool(int(sys.argv[5]))

print("HESE Pulses: {}, Use FADC {}, Use ATWD {}, Beacon FADC {}, Noise Cut {}, Causal Qtot {}".format(hese_pulses, use_fadc, use_atwd, beacon_fadc, noise_cut, causal_qtot))
print("----------------------------")

plot_causal_homog_comparison=True
if plot_causal_homog_comparison:

	data_homoq = h5py.File('comparison_overlap_HESE_SplitInIcePulses_FADC_True_ATWD_True_FADCBeacon_True_NoiseCut_True_CausalQtot_False.hdf5', 'r')
	data_causq = h5py.File('comparison_overlap_HESE_SplitInIcePulses_FADC_True_ATWD_True_FADCBeacon_True_NoiseCut_True_CausalQtot_True.hdf5', 'r')


	'''
		This bit deals with making plots of what happened to the "overlap" DOMS
		That is, what changed about the DOMs that are in both EHE and HESE
	'''
	plot_overlap = True
	if plot_overlap:
		strings = np.asarray(data_homoq['data_overlap']['strings'])
		doms = np.asarray(data_homoq['data_overlap']['doms'])
		homoq_char = np.asarray(data_homoq['data_overlap']['hese_overlap'])
		causq_char = np.asarray(data_causq['data_overlap']['hese_overlap'])

		print("homog qtot overlap {}".format(np.sum(homoq_char)))
		print("causal qtot overlap {}".format(np.sum(causq_char)))

		total_correction = 0
		corrected = []
		for s, o, h, c in zip(strings, doms, homoq_char, causq_char):
			if h-c > 5:
				print("{}, {}: Homog qtot {:6.2f}, Causal qtot {:6.2f}, diff {:6.2f}".format(s, o, h, c, h-c))
				total_correction+=h-c
				corrected.append(h-c)
		print("Total just found in >0 diff {} in {} corrected DOMs".format(total_correction, len(corrected)))

		fig, axs = plt.subplots(1,1,figsize=(5,5))
		axs.plot(causq_char, homoq_char, 'o', alpha=0.5)
		axs.plot([0,1000],[0,1000],'--', label='1-1 line')
		axs.set_xlabel('Causal Qtot [pe]')
		axs.set_ylabel('Homogenized Qtot [pe]')
		axs.set_xlim([-10, 250])
		axs.set_ylim([-10, 250])
		axs.set_aspect('equal')
		# axs.set_xlim([-1,10])
		# axs.set_ylim([-1,10])
		plt.tight_layout()
		fig.savefig('homog_causal_qtot_comparison_overlap.png', dpi=300)

		# and a histogram
		fig, axs = plt.subplots(1,2,figsize=(10,5))
		# bins = np.linspace(-0.15,0.15,100)
		bins=50
		axs[0].hist(causq_char, bins=bins, alpha=0.5, label='Causal')#, histtype='step')
		axs[0].hist(homoq_char, bins=bins, alpha=0.5, label='Homogenized')#, histtype='step')
		axs[0].set_xlabel('Charge [pe]')
		axs[0].set_ylabel('Number of DOMs')
		axs[0].set_yscale('log')
		axs[0].legend()

		correction = homoq_char - causq_char
		mask = correction>0
		# bins = np.linspace(-1,25,52)
		bins = 50
		stuff = axs[1].hist((homoq_char - causq_char)[mask], color='C2', bins=bins, alpha=0.5)
		axs[1].set_xlabel('Homogenized - Causal [pe]')
		axs[1].set_ylabel('Number of DOMs')
		# axs[1].set_yscale('log')
		# axs[1].plot([2,2], [1,20], 'k--')
		print("total q difference {}".format(np.sum(homoq_char - causq_char)))

		plt.tight_layout()
		fig.savefig('homog_causal_qtot_comparison_hist_overlap.png', dpi=300)
		del fig, axs


	'''
		This bit deals with making plots of what happened to the "excluded" DOMS
		That is, what changed about the DOMs that are *only* in HESE
	'''

	plot_exclusion = False
	if plot_exclusion:
		homoq_keys = np.asarray(data_homoq['data_exclusion']['missing_in_ehe_keys'])
		homoq_char = np.asarray(data_homoq['data_exclusion']['missing_in_ehe'])
		causq_keys = np.asarray(data_causq['data_exclusion']['missing_in_ehe_keys'])
		causq_char = np.asarray(data_causq['data_exclusion']['missing_in_ehe'])
		mask_hq = homoq_char>0
		mask_cq = causq_char>0
		print('Number of DOMs exclusively in HESE in Hqtot calculation: {}'.format(len(homoq_char[mask_hq])))
		print('Number of DOMs exclusively in HESE in Cqtot calculation: {}'.format(len(causq_char[mask_cq])))
		
		new_doms = []
		old_doms_difference = []
		for key, hq, cq in zip(homoq_keys, homoq_char, causq_char):
			if not cq>0:
				print("{}: Hqtot {:.2f}, Cqtot {:.2f}".format(key, hq, cq))
				new_doms.append(hq)
			else:
				old_doms_difference.append(hq-cq)

		print("{} pe of charge is found in {} new DOMs".format(np.sum(new_doms), len(new_doms)))
		print("{} pe of charge is found in adjustments to {} old DOMs".format(np.sum(old_doms_difference), len(old_doms_difference)))

		# hist
		fig, axs = plt.subplots(1,1,figsize=(5,5))
		bins=np.linspace(0,20,20)
		axs.hist(new_doms, bins=bins, alpha=0.5, label='{:.2f} pe in {} New DOMs'.format(np.sum(new_doms), len(new_doms)))#, histtype='step')
		axs.hist(old_doms_difference, bins=bins, alpha=0.5, label='{:.2f} pe in {} Old DOMs'.format(np.sum(old_doms_difference), len(old_doms_difference)))#, histtype='step')
		axs.set_xlabel(r'$Q_{homogenized} - Q_{causal}$ [pe]')
		axs.set_ylabel('Number of DOMs')
		axs.set_yscale('log')
		axs.legend()
		plt.tight_layout()
		fig.savefig('homog_causal_qtot_changes_hist.png', dpi=300)
		del fig, axs

		fig, axs = plt.subplots(1,1,figsize=(5,5))
		bins = np.logspace(np.log10(0.1), np.log10(2e1), 50)
		axs.hist(homoq_char, bins=bins, alpha=1.0, color='C0', label='Homogenized: {:.2f} pe in {} DOMs'.format(np.sum(homoq_char), len(homoq_char[mask_hq])), histtype='step', linewidth=2)
		axs.hist(causq_char, bins=bins, alpha=0.5, color='C1', label='Causal: {:.2f} pe in {} DOMs'.format(np.sum(causq_char), len(causq_char[mask_cq])))#, histtype='step', linewidth=2)
		axs.set_title('Charge Missing in EHE Estimate')
		axs.set_xlabel(r'NPE')
		axs.set_ylabel('Number of DOMs')
		axs.set_ylim([0.9,3e1])
		axs.set_yscale('log')
		axs.set_xscale('log')
		axs.legend()
		plt.tight_layout()
		fig.savefig('homog_causal_qtot_missing_in_ehe.png', dpi=300)
		del fig, axs	

	do_old = False # just want to be able to collapse and ignore this (it's not as useful as I was hoping)
	if do_old:

		baselines_in = h5py.File('comparison_causal_homog_qtot_SplitInIcePulses.hdf5', 'r')
		data = baselines_in['charge']
		strings = data['strings']
		doms = data['doms']

		causal_pe = np.asarray(data['causal_pe'])
		homog_pe = np.asarray(data['homog_pe'])

		fig, axs = plt.subplots(1,1,figsize=(5,5))
		axs.plot(causal_pe, homog_pe, 'o', alpha=0.5)
		axs.plot([0,1000],[0,1000],'--', label='1-1 line')
		axs.set_xlabel('Causal Qtot [pe]')
		axs.set_ylabel('Homogenized Qtot [pe]')
		axs.set_xlim([-10, 250])
		axs.set_ylim([-10, 250])
		axs.set_aspect('equal')
		# axs.set_xlim([-1,10])
		# axs.set_ylim([-1,10])
		plt.tight_layout()
		fig.savefig('homog_causal_qtot_comparison.png', dpi=300)

		# axs.set_xlim([0.1, 1000])
		# axs.set_ylim([0.1, 1000])
		# axs.set_xscale('log')
		# axs.set_yscale('log')
		# axs.set_aspect('equal')
		# plt.tight_layout()
		# fig.savefig('homog_causal_qtot_comparison_log.png', dpi=300)

		del fig, axs

		# fig, axs = plt.subplots(1,1,figsize=(5,5))
		# my_map = plt.cm.plasma
		# counts, xedges, yedges, im = axs.hist2d(causal_pe, homog_pe, 
		# 	bins=50,
		# 	range = [[-1,10],[-1,10]],
		# 	cmap=my_map,
		# 	# norm=colors.LogNorm(),
		# 	cmin=1)
		# cbar = plt.colorbar(im, ax=axs)
		# cbar.set_label('Number of DOMs')
		# axs.set_xlabel(r'Causal Qtot [pe]')
		# axs.set_ylabel(r'Homogenized Qtot [pe]')
		# axs.set_aspect('equal')	
		# plt.tight_layout()
		# fig.savefig('homog_causal_qtot_comparison_2dhist.png', dpi=300)


		# and a histogram
		fig, axs = plt.subplots(1,2,figsize=(10,5))
		# bins = np.linspace(-0.15,0.15,100)
		bins=50
		axs[0].hist(causal_pe, bins=bins, alpha=0.5, label='Causal')#, histtype='step')
		axs[0].hist(homog_pe, bins=bins, alpha=0.5, label='Homogenized')#, histtype='step')
		axs[0].set_xlabel('Charge [pe]')
		axs[0].set_ylabel('Number of DOMs')
		axs[0].set_yscale('log')
		axs[0].legend()

		bins = np.linspace(-1,25,52)
		axs[1].hist(homog_pe - causal_pe, color='C2', bins=bins, alpha=0.5)
		axs[1].set_xlabel('Homogenized - Causal [pe]')
		axs[1].set_ylabel('Number of DOMs')
		axs[1].set_yscale('log')
		print("total q difference {}".format(np.sum(homog_pe - causal_pe)))

		plt.tight_layout()
		fig.savefig('homog_causal_qtot_hist.png', dpi=300)
		del fig, axs


plot_baseline_stuff = False
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

	file_in = h5py.File('comparison_overlap_HESE_{}_FADC_{}_ATWD_{}_FADCBeacon_{}_NoiseCut_{}_CausalQtot_{}.hdf5'.format(hese_pulses, use_fadc, use_atwd, beacon_fadc, noise_cut, causal_qtot), 'r')
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
		this_ratio = this_hese/this_ehe
		# print("({}, {}): EHE {:.2f}, HESE {:.2f}, rel dif {:.2f}".format(string, dom, this_ehe, this_hese, rel_diff))	
		if this_ratio > 1.15:
			# print("({}, {}): EHE {:.2f}, HESE {:.2f}, dif {:.2f}, ratio {:.2f}".format(string, dom, this_ehe, this_hese, this_hese-this_ehe, this_hese/this_ehe))
			big_ehe.append(this_ehe)
			big_hese.append(this_hese)
		diff.append(rel_diff)
		ratio.append(this_ehe/this_hese)
	diff = np.asarray(diff)
	big_ehe = np.asarray(big_ehe)
	big_hese = np.asarray(big_hese)
	ratio = np.asarray(ratio)


	# make a plot of EHE Q vs HESE Q
	fig, axs = plt.subplots(1,1,figsize=(5,5))
	print("Number of mutual DOMs is {}".format(len(ehe_npe)))
	axs.plot(ehe_npe, hese_npe,'o', alpha=0.5, label='N={}'.format(len(ehe_npe)))
	axs.set_xlabel('EHE NPE')
	axs.set_ylabel('HESE NPE ({})'.format(hese_pulses))
	axs.plot([0,300],[0,300],'--', label='1-1 line')
	# axs.set_title('FADC = {}, ATWD = {}, FADC Beacon Baselines {}, Noise Cut'.format(use_fadc, use_atwd, beacon_fadc, noise_cut))
	axs.set_title('Use FADC = {}, Use ATWD = {}, \nUse FADC Beacon Baselines={}, \n Use FADC Noise Cut={}, Causal Qtot={}'.format(use_fadc, use_atwd, beacon_fadc, noise_cut, causal_qtot), size=10)
	axs.legend()
	# axs.plot(big_ehe, big_hese, 'o', alpha=0.5, color='red')

	print("Total charge in overlapping DOMs in HESE is {:.2f}".format(np.sum(hese_npe)))
	print("Total charge in overlapping DOMs in EHE  is {:.2f}".format(np.sum(ehe_npe)))


	do_log = False
	if not do_log:
		axs.set_xlim([-10, 300])
		axs.set_ylim([-10, 300])
		axs.set_aspect('equal')
		axs.legend()
		plt.tight_layout()
		fig.savefig('ehenpe_vs_hesenpe_HESE_{}_FADC_{}_ATWD_{}_FADCBeacon_{}_NoiseCut_{}_CausalQtot_{}.png'.format(hese_pulses, use_fadc, use_atwd, beacon_fadc, noise_cut, causal_qtot), dpi=300)
		plt.close(fig)
		del fig, axs
	else:
		axs.set_yscale('log')
		axs.set_xscale('log')
		axs.set_xlim([0.1, 300])
		axs.set_ylim([0.1, 300])
		axs.set_aspect('equal')
		plt.tight_layout()
		fig.savefig('ehenpe_vs_hesenpe_HESE_{}_FADC_{}_ATWD_{}_FADCBeacon_{}_NoiseCut_{}_CausalQtot_{}_log.png'.format(hese_pulses, use_fadc, use_atwd, beacon_fadc, noise_cut, causal_qtot), dpi=300)
		plt.close(fig)
		del fig, axs

	bins = np.linspace(0,25,26)
	fig, axs = plt.subplots(1,1,figsize=(5,5))
	axs.hist(ratio, bins=bins, alpha=0.5, label='N={}'.format(len(ratio)))#, histtype='step')
	# axs.set_yscale('log')
	axs.set_xlabel('EHE/HESE  [HESE={}]'.format(hese_pulses))
	axs.set_ylabel('Number of DOMs')
	# axs.set_xlim([-1.0, 1.0])
	axs.set_title('Use FADC = {}, Use ATWD = {}, \nUse FADC Beacon Baselines={}, \n Use FADC Noise Cut={}, Causal Qtot={}'.format(use_fadc, use_atwd, beacon_fadc, noise_cut, causal_qtot), size=10)
	# axs.set_ylim([0,60])
	# axs.plot([0,250],[0,250],'--', label='1-1 line')
	axs.set_yscale('log')
	axs.legend()
	plt.tight_layout()
	# fig.savefig('ehe_over_hese_{}_FADC_{}_ATWD_{}_forcezero_{}.png'.format(hese_pulses, use_fadc, use_atwd, force_zero), dpi=300)
	fig.savefig('ehenpe_over_hese_HESE_{}_FADC_{}_ATWD_{}_FADCBeacon_{}_NoiseCut_{}_CausalQtot_{}.png'.format(hese_pulses, use_fadc, use_atwd, beacon_fadc, noise_cut, causal_qtot), dpi=300)
	del fig, axs

	##############################
	##############################

	fig, axs = plt.subplots(1,1,figsize=(5,5))
	axs.hist(diff, bins=50, alpha=0.5, label='N={}'.format(len(diff)))#, histtype='step')
	# axs.set_yscale('log')
	axs.set_xlabel('(HESE - EHE)/EHE  [HESE={}]'.format(hese_pulses))
	axs.set_ylabel('Number of DOMs')
	# axs.set_xlim([-1.0, 1.0])
	axs.set_yscale('log')
	axs.set_title('Use FADC = {}, Use ATWD = {}, \nUse FADC Beacon Baselines={}, \n Use FADC Noise Cut={}, Causal Qtot={}'.format(use_fadc, use_atwd, beacon_fadc, noise_cut, causal_qtot), size=10)
	# axs.set_ylim([0,60])
	# axs.plot([0,250],[0,250],'--', label='1-1 line')
	axs.legend()
	plt.tight_layout()
	# fig.savefig('ehe_vs_hese_rel_diff_{}_FADC_{}_ATWD_{}_forcezero_{}.png'.format(hese_pulses, use_fadc, use_atwd, force_zero), dpi=300)
	fig.savefig('ehe_vs_hese_reldiff_HESE_{}_FADC_{}_ATWD_{}_FADCBeacon_{}_NoiseCut_{}_CausalQtot_{}.png'.format(hese_pulses, use_fadc, use_atwd, beacon_fadc, noise_cut, causal_qtot), dpi=300)

	##############################
	##############################	

	data_exclusion = file_in['data_exclusion']
	missing_in_hese = np.asarray(data_exclusion['missing_in_hese'])
	missing_in_ehe = np.asarray(data_exclusion['missing_in_ehe'])
	mask = missing_in_ehe > 0
	print("Total missing in ehe {} in {} doms".format(np.sum(missing_in_ehe), len(missing_in_ehe[mask])))

	fig, axs = plt.subplots(1,1,figsize=(5,5))
	# bins = range(-10,1000+10,10)
	bins = np.logspace(np.log10(0.1), np.log10(1e2), 50)
	# stuff = axs.hist(missing_in_hese, bins=bins, alpha=0.5, label='Missing in HESE ({:.2f} PE)'.format(np.sum(missing_in_hese)))#, histtype='step')
	stuff2 = axs.hist(missing_in_ehe, bins=bins, alpha=0.5, label='Missing in EHE ({:.2f} PE in {} DOMs)'.format(np.sum(missing_in_ehe), len(missing_in_ehe[mask])))#, histtype='step')
	axs.set_title('Use FADC = {}, Use ATWD = {}, \nUse FADC Beacon Baselines={}, \n Use FADC Noise Cut={}, Causal Qtot={}'.format(use_fadc, use_atwd, beacon_fadc, noise_cut, causal_qtot), size=10)
	axs.set_xlabel('NPE')
	axs.set_ylabel('Number of DOMs')
	axs.set_yscale('log')
	axs.set_xscale('log')
	axs.set_ylim([0.1,100])
	axs.legend()
	plt.tight_layout()
	# fig.savefig('ehe_vs_hese_missing_{}_FADC_{}_ATWD_{}_forcezero_{}.png'.format(hese_pulses, use_fadc, use_atwd, force_zero), dpi=300)
	fig.savefig('ehe_vs_hese_missing_reldiff_HESE_{}_FADC_{}_ATWD_{}_FADCBeacon_{}_NoiseCut_{}_CausalQtot_{}.png'.format(hese_pulses, use_fadc, use_atwd, beacon_fadc, noise_cut, causal_qtot), dpi=300)

