import numpy as np
import h5py as h5py
import matplotlib.pyplot as plt

file_in = h5py.File('comparison_overlap.hdf5', 'r')
data_overlap = file_in['data_overlap']
hese_npe = np.asarray(data_overlap['hese_overlap'])
ehe_npe = np.asarray(data_overlap['ehe_overlap'])
strings = np.asarray(data_overlap['strings'])
doms = np.asarray(data_overlap['doms'])

# make a plot of EHE Q vs HESE Q
fig, axs = plt.subplots(1,1,figsize=(5,5))
axs.plot(ehe_npe, hese_npe,'o', alpha=0.5, label='data')
axs.set_xlabel('EHE NPE')
axs.set_ylabel('HESE NPE')
axs.plot([0,300],[0,300],'--', label='1-1 line')

axs.set_xlim([-10, 300])
axs.set_ylim([-10, 300])
axs.legend()
axs.set_aspect('equal')
plt.tight_layout()
fig.savefig('ehenpe_vs_hesenpe.png', dpi=300)

axs.set_yscale('log')
axs.set_xscale('log')
axs.set_xlim([0.1, 300])
axs.set_ylim([0.1, 300])
axs.set_aspect('equal')
plt.tight_layout()
# fig.savefig('ehenpe_vs_hesenpe_log.png', dpi=300)
plt.close(fig)
del fig, axs


# plot their differences
diff = []
for this_ehe, this_hese in zip(ehe_npe, hese_npe):
	rel_diff = (this_hese - this_ehe)/this_ehe
	diff.append(rel_diff)
diff = np.asarray(diff)

fig, axs = plt.subplots(1,1,figsize=(5,5))
axs.hist(diff, bins=50, alpha=0.5)#, histtype='step')
# axs.set_yscale('log')
axs.set_xlabel('(HESE - EHE)/EHE')
axs.set_ylabel('Number of DOMs')
axs.set_xlim([-1.0, 1.0])
# axs.plot([0,250],[0,250],'--', label='1-1 line')
# axs.legend()
plt.tight_layout()
fig.savefig('ehe_vs_hese_rel_diff.png', dpi=300)

##############################
##############################

data_exclusion = file_in['data_exclusion']
missing_in_hese = np.asarray(data_exclusion['missing_in_hese'])
missing_in_ehe = np.asarray(data_exclusion['missing_in_ehe'])

fig, axs = plt.subplots(1,1,figsize=(5,5))
# bins = range(-10,1000+10,10)
bins = np.logspace(np.log10(0.1), np.log10(1e3), 50)
axs.hist(missing_in_hese, bins=bins, alpha=0.5, label='Missing in HESE')#, histtype='step')
axs.hist(missing_in_ehe, bins=bins, alpha=0.5, label='Missing in EHE')#, histtype='step')
axs.set_xlabel('NPE')
axs.set_ylabel('Number of DOMs')
# axs.set_yscale('log')
axs.set_xscale('log')
axs.legend()
plt.tight_layout()
fig.savefig('ehe_vs_hese_missing.png', dpi=300)
