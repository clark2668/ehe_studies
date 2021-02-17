import numpy as np
import h5py as h5py
import matplotlib.pyplot as plt

use_fadc=False
hese_pulses = 'SplitInIcePulses'

print("HESE Pulses: {}, FADC: {}".format(hese_pulses, use_fadc))
print("----------------------------")

file_in = h5py.File('comparison_overlap_{}_FADC_{}.hdf5'.format(hese_pulses, use_fadc), 'r')
data_overlap = file_in['data_overlap']
hese_npe = np.asarray(data_overlap['hese_overlap'])
ehe_npe = np.asarray(data_overlap['ehe_overlap'])
strings = np.asarray(data_overlap['strings'])
doms = np.asarray(data_overlap['doms'])


# plot their differences
diff = []
big_ehe = []
big_hese = []
for this_ehe, this_hese, string, dom in zip(ehe_npe, hese_npe, strings, doms):
	rel_diff = (this_hese - this_ehe)/this_ehe
	if rel_diff > .3:
		print("({}, {}): this_ehe {:.2f}, this_hese {:.2f}, rel dif {:.2f}".format(string, dom, this_ehe, this_hese, rel_diff))
		big_ehe.append(this_ehe)
		big_hese.append(this_hese)
	diff.append(rel_diff)
diff = np.asarray(diff)
big_ehe = np.asarray(big_ehe)
big_hese = np.asarray(big_hese)

# make a plot of EHE Q vs HESE Q
fig, axs = plt.subplots(1,1,figsize=(5,5))
print("Number of mutual DOMs is {}".format(len(ehe_npe)))
axs.plot(ehe_npe, hese_npe,'o', alpha=0.5, label='N={}'.format(len(ehe_npe)))
axs.set_xlabel('EHE NPE')
axs.set_ylabel('HESE NPE ({})'.format(hese_pulses))
axs.plot([0,300],[0,300],'--', label='1-1 line')
axs.set_title('Include FADC = {}'.format(use_fadc))
axs.legend()
# axs.plot(big_ehe, big_hese, 'o', alpha=0.5, color='red')

axs.set_xlim([-10, 300])
axs.set_ylim([-10, 300])
axs.legend()
axs.set_aspect('equal')
plt.tight_layout()
fig.savefig('ehenpe_vs_hesenpe_{}_FADC_{}.png'.format(hese_pulses, use_fadc), dpi=300)

# axs.set_yscale('log')
# axs.set_xscale('log')
# axs.set_xlim([0.1, 300])
# axs.set_ylim([0.1, 300])
# axs.set_aspect('equal')
# plt.tight_layout()
# fig.savefig('ehenpe_vs_hesenpe_{}_FADC_{}_log.png'.format(hese_pulses, use_fadc), dpi=300)
plt.close(fig)
del fig, axs



fig, axs = plt.subplots(1,1,figsize=(5,5))
axs.hist(diff, bins=50, alpha=0.5, label='N={}'.format(len(diff)))#, histtype='step')
# axs.set_yscale('log')
axs.set_xlabel('(HESE - EHE)/EHE  [HESE={}]'.format(hese_pulses))
axs.set_ylabel('Number of DOMs')
axs.set_xlim([-1.0, 1.0])
axs.set_title('Include FADC = {}'.format(use_fadc))
axs.set_ylim([0,60])
# axs.plot([0,250],[0,250],'--', label='1-1 line')
axs.legend()
plt.tight_layout()
fig.savefig('ehe_vs_hese_rel_diff_{}_FADC_{}.png'.format(hese_pulses, use_fadc), dpi=300)

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
axs.set_title('Include FADC = {}'.format(use_fadc))
axs.set_xlabel('NPE')
axs.set_ylabel('Number of DOMs')
# axs.set_yscale('log')
axs.set_xscale('log')
axs.legend()
plt.tight_layout()
fig.savefig('ehe_vs_hese_missing_{}_FADC_{}.png'.format(hese_pulses, use_fadc), dpi=300)
