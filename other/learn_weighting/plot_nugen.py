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
weights_atmo = np.asarray([])
weights_astro = np.asarray([])
weights_cosmo = np.asarray([])
weights = np.asarray([])
npe = np.asarray([])
zenith = np.asarray([])

sim_weights = None
weighter = None

data = np.genfromtxt("ahlers_100.txt",delimiter='\t',skip_header=8,names=['e','f', 'no1', 'no2'])
# spline in GeV and log10(flux) to make things easier
# interpolator = splrep(np.log10(data['e'])-9, np.log10(data['f']))
interpolator = interp1d(np.log10(data['e'])-9, np.log10(data['f']),
	bounds_error=False, fill_value=np.nan)

def ahlers(energy):
	# energy in GeV
	this_e = np.log10(energy)
	# return 10**splev(this_e, interpolator)
	return (10**interpolator(this_e))/np.power(energy, 2.)

def power_law(energy):
	# https://pos.sissa.it/301/1005/pdf
	return 1.01e-18 * np.power(energy/1e5,-2.19)

cr_model = nuflux.makeFlux('H3a_SIBYLL23C')

for file in files:
	file_in = pd.HDFStore(file, 'r')
	print(file)
	if args.type == 'nugen':
		loc_weighter = simweights.NuGenWeighter(file_in, nfiles=1000)
		flux = power_law
	elif args.type == 'corsika':
		loc_weighter = simweights.CorsikaWeighter(file_in, nfiles=1000)
		flux = simweights.GaisserH4a()
	else:
		raise AttributeError('Type {} is not supported'.format(args.type))
	
	local_primary_energy = loc_weighter.get_column('PolyplopiaPrimary', 'energy')
	local_weights = loc_weighter.get_weights(flux)
	local_npe = loc_weighter.get_column('EHEPortiaEventSummarySRT','bestNPEbtw')
	local_zenith = loc_weighter.get_column('EHEOpheliaParticleSRT_ImpLF', 'zenith')

	local_weights_astro = loc_weighter.get_weights(power_law)
	local_weights_atmo = loc_weighter.get_weights(cr_model)
	local_weights_cosmo = loc_weighter.get_weights(ahlers)
	weights_astro = np.concatenate((weights_astro, local_weights_astro))
	weights_atmo = np.concatenate((weights_atmo, local_weights_atmo))
	weights_cosmo = np.concatenate((weights_cosmo, local_weights_cosmo))

	primary_energy = np.concatenate((primary_energy, local_primary_energy))
	weights = np.concatenate((weights, local_weights))
	npe = np.concatenate((npe, local_npe))
	zenith = np.concatenate((zenith, local_zenith))

	file_in.close()

sizer=15

# # 2D hist
# fig = plt.figure(figsize=(9,7))
# ax = fig.add_subplot(111)
# counts, xedges, yedges, im = ax.hist2d(np.cos(zenith), np.log10(npe), 
# 		bins=100,
# 		# range = [[plot_min,plot_max],[plot_min,plot_max]],
# 		cmap=plt.cm.plasma,
# 		norm=colors.LogNorm(),
# 		# cmin=1
# 		)
# cbar = plt.colorbar(im, ax=ax)
# # ax.plot([plot_min,plot_max],[plot_min,plot_max], 'k:')
# cbar.set_label('Weighted Number of Events', fontsize=sizer)
# ax.set_ylabel(r'log$_{10}$(Portia NPE)', fontsize=sizer)
# ax.set_xlabel(r'cos$\theta$', fontsize=sizer)
# ax.tick_params(labelsize=sizer)
# # plt.tight_layout()
# fig.savefig('npe_vs_coszen.png', edgecolor='none', bbox_inches='tight')
# del fig, ax
# hist of rate vs true energy
fig = plt.figure(figsize=(5,6))
ax = fig.add_subplot(111)
bins = np.geomspace(1e5, 1e9, 50)
ax.hist(primary_energy, bins=bins, 
	histtype='step', alpha=0.75, label=r'Unweighted (NuMu, 100 TeV - 500 PeV, E$^{-1}$)')
ax.hist(primary_energy, weights=weights_astro, bins=bins, 
	histtype='step', alpha=0.75, label=r'Astro, E$^{-2.19}$, PoS(ICRC2017)1005')
ax.hist(primary_energy, weights=weights_atmo, bins=bins, 
	histtype='step', alpha=0.75, label='Atmo, Gaisser H3a w/ SIBYLL23C')
nan_mask = ~np.isnan(weights_cosmo) # mask the out of range bits of the ahlers flux
ax.hist(primary_energy[nan_mask], weights=weights_cosmo[nan_mask], bins=bins, 
	histtype='step', alpha=0.75, label='Cosmo, Ahlers 100% proton')
ax.set_ylabel('Rate [Hz]')
ax.set_xlabel('True Neutrino Energy [GeV]')
ax.set_yscale('log')
ax.set_xscale('log')
ax.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
                mode="expand", borderaxespad=0)#, ncol=3)
plt.tight_layout()
fig.savefig('weighted_energy_distribution.png', dpi=300)
del fig, ax
