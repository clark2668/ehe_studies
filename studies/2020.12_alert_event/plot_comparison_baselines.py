import numpy as np
import h5py as h5py
import matplotlib.pyplot as plt

which_atwd = 'atwd_0_1'

file_in = h5py.File('compare_baselines.hdf5', 'r')
data = file_in['baselines']
atwd_0_0 = np.asarray(data[which_atwd])
portia = np.asarray(data['portia'])
strings = np.asarray(data['strings'])
doms = np.asarray(data['doms'])

fig, axs = plt.subplots(1,1,figsize=(5,5))
axs.plot(portia, strings,'o', alpha=0.5)
axs.set_xlabel(f'Portia Baseline')
axs.set_ylabel(f'Beacon Baseline {which_atwd}')
# axs.plot([0,300],[0,300],'--', label='1-1 line')
# axs.set_xlim([-10, 300])
# axs.set_ylim([-10, 300])
# axs.legend()
# axs.set_aspect('equal')
plt.tight_layout()
fig.savefig(f'portiabaseline_vs_beaconbaseline_{which_atwd}.png', dpi=300)
