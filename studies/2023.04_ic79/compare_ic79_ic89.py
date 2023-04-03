import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import style
import tables
import pandas as pd
import copy
import yaml
import numpy as np


livetime = 365*24*60*60

cfg_file = 'config.yaml'
cfg_file = yaml.safe_load(open(cfg_file))


def harvest_values(cfg_file, dataset):
    return_dict = {}
    the_f = tables.open_file(cfg_file['burn_sample'][dataset]['file'])
    charge = the_f.get_node(f"/{cfg_file['variables']['charge']['variable']}").col(f"{cfg_file['variables']['charge']['value']}")
    return_dict['charge'] = copy.deepcopy(charge)
    return_dict['livetime'] = cfg_file['burn_sample'][dataset]['livetime']
    the_f.close()
    qmask = return_dict['charge'] > 1E3
    return_dict['qmask'] = qmask
    return return_dict

ic79_dict = harvest_values(cfg_file, 'IC79-2010-pass2')
ic86_dict = harvest_values(cfg_file, 'IC86-2011-pass2')

fig, (ax, axr) = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [3, 1]})

q_bins = np.logspace(4,6,10)
def get_bin_centers(bins):
    return (bins[1:] + bins[:-1]) * 0.5
qbin_centers = var_bin_centers = get_bin_centers(q_bins)

ic79_sum, b = np.histogram(ic79_dict['charge'][ic79_dict['qmask']],bins=q_bins)
ic86_sum, b = np.histogram(ic86_dict['charge'][ic86_dict['qmask']],bins=q_bins)

# rescale to 1 year
ic79_sum, b, p = ax.hist(
    x=qbin_centers,
    weights=ic79_sum/ic79_dict['livetime']*livetime,
    bins=q_bins,
    label='IC79',
    histtype='step', lw=2
)

ic86_sum, b, p = ax.hist(
    x=qbin_centers,
    weights=ic86_sum/ic86_dict['livetime']*livetime,
    bins=q_bins,
    label='IC86',
    histtype='step', lw=2
)

axr.plot(
    qbin_centers,
    ic79_sum/ic86_sum,
    'o-'
)

ax.legend()
ax.set_yscale('log')
ax.set_xscale('log')
# ax.set_ylim([1E-5, 1E4])
ax.set_ylabel('Events / {:.2f} days (rescaled)'.format(livetime/(60*60*24)))
axr.set_ylabel('IC79/IC86')
axr.set_xlabel(r'Charge')
# axr.set_ylim([0.5,1.5])
axr.grid()
axr.axhline(y=1,linestyle='--', color='red')
fig.tight_layout()
fig.savefig('plots/charge.png')
del fig, ax, axr
