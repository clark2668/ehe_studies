import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import argparse
from scipy.interpolate import splrep, splev, interp1d

import simweights
import nuflux


parser = argparse.ArgumentParser()
parser.add_argument("-f", type=str, nargs='+',
	dest="input_files", required=True,
	help="paths to input files (absolute or relative)")
parser.add_argument("-t", type=str, 
	dest="type", required=True,
	help="weighter type (nugen or corsika)")
args = parser.parse_args()
files = args.input_files

primary_energy = np.asarray([])
weights = np.asarray([])
npe = np.asarray([])
zenith = np.asarray([])

sim_weights = None
weighter = None

for file in files:
	file_in = pd.HDFStore(file, 'r')
	print(file)
	if args.type == 'nugen':
		loc_weighter = simweights.NuGenWeighter(file_in, nfiles=1000)
		flux = power_law
	elif args.type == 'corsika':
		loc_weighter = simweights.CorsikaWeighter(file_in, nfiles=int(len(files)))
		flux = simweights.GaisserH4a()
	else:
		raise AttributeError('Type {} is not supported'.format(args.type))
	
	local_primary_energy = loc_weighter.get_column('PolyplopiaPrimary', 'energy')
	local_weights = loc_weighter.get_weights(flux)
	local_npe = loc_weighter.get_column('EHEPortiaEventSummarySRT','bestNPEbtw')
	local_zenith = loc_weighter.get_column('EHEOpheliaParticleSRT_ImpLF', 'zenith')

	primary_energy = np.concatenate((primary_energy, local_primary_energy))
	weights = np.concatenate((weights, local_weights))
	npe = np.concatenate((npe, local_npe))
	zenith = np.concatenate((zenith, local_zenith))

	file_in.close()

sizer=15

# 2D hist
fig = plt.figure(figsize=(18,7))
ax = fig.add_subplot(121)
bins = [np.linspace(-1,1, 50), np.linspace(2, 8, 60)]
counts, xedges, yedges, im = ax.hist2d(np.cos(zenith), np.log10(npe), 
		bins=bins,
		cmap=plt.cm.plasma,
		norm=colors.LogNorm(),
		# cmin=1
		)
cbar = plt.colorbar(im, ax=ax)
# ax.plot([plot_min,plot_max],[plot_min,plot_max], 'k:')
cbar.set_label('Number of Thrown Events', fontsize=sizer)
ax.set_ylabel(r'log$_{10}$(Portia NPE)', fontsize=sizer)
ax.set_xlabel(r'cos$\theta$', fontsize=sizer)
ax.tick_params(labelsize=sizer)

ax2 = fig.add_subplot(122)
counts, xedges, yedges, im = ax2.hist2d(np.cos(zenith), np.log10(npe), 
		bins=bins,
		cmap=plt.cm.plasma,
		norm=colors.LogNorm(),
		weights=weights,
		# cmin=1E-7
		)
cbar2 = plt.colorbar(im, ax=ax2)
cbar2.set_label('Rate [Hz]', fontsize=sizer)
ax2.set_ylabel(r'log$_{10}$(Portia NPE)', fontsize=sizer)
ax2.set_xlabel(r'cos$\theta$', fontsize=sizer)
ax2.tick_params(labelsize=sizer)

plt.tight_layout()
fig.savefig('npe_vs_coszen_corsika.png', edgecolor='none', bbox_inches='tight')
del fig, ax

# hist of rate vs true energy
fig = plt.figure(figsize=(5,6))
ax = fig.add_subplot(111)
bins = np.geomspace(1e5, 1e10, 50)
ax.hist(primary_energy, bins=bins, 
	histtype='step', alpha=0.75, label=r'Unweighted (Corsika, 1 PeV - 10 EeV, E$^{-2}$)')
ax.hist(primary_energy, weights=weights, bins=bins, 
	histtype='step', alpha=0.75, label=r'Weighted to Gaisser H4a')
ax.set_ylabel('Rate [Hz]')
ax.set_xlabel('Energy [GeV]')
ax.set_yscale('log')
ax.set_xscale('log')
ax.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
                mode="expand", borderaxespad=0)#, ncol=3)
plt.tight_layout()
fig.savefig('weighted_energy_distribution_corsika.png', dpi=300)
del fig, ax
