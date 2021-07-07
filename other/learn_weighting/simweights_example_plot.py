import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse

import simweights
import nuflux


parser = argparse.ArgumentParser()
parser.add_argument("-f", type=str, nargs='+',
	dest="input_files", required=True,
	help="paths to input files (absolute or relative)")
args = parser.parse_args()
files = args.input_files

primary_energy = np.asarray([])
weights = np.asarray([])

# weighter = None

def power_law(energy):
	# https://pos.sissa.it/301/1005/pdf
	return 1.01e-18 * np.power(energy/1e5,-2.19)

cr_model = nuflux.makeFlux('H3a_SIBYLL23C')

for file in files:
	file_in = pd.HDFStore(file, 'r')
	loc_weighter = simweights.NuGenWeighter(file_in, nfiles=10)

	# if weighter is None:
	# 	weighter = loc_weighter
	# else:
	# 	weighter += loc_weighter
	
	local_primary_energy = loc_weighter.get_column('PolyplopiaPrimary', 'energy')
	# local_weights = loc_weighter.get_weights(power_law)
	local_weights = loc_weighter.get_weights(cr_model)
	
	primary_energy = np.concatenate((primary_energy, local_primary_energy))
	weights = np.concatenate((weights, local_weights))

	file_in.close()

print("Total rate: {} Hz".format(np.sum(weights)))

fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(111)
bins = np.geomspace(1e2, 1e8, 50)
# ax.hist(primary_energy, bins=bins, alpha=0.5, label='Unweighted')
ax.hist(primary_energy, weights=weights, bins=bins, alpha=0.5, label='Weighted')
ax.set_ylabel('Rate [Hz]')
ax.set_xlabel('Energy [GeV]')
ax.set_yscale('log')
ax.set_xscale('log')
ax.legend()
plt.tight_layout()
fig.savefig('weighted_energy_distribution.png', dpi=300)
del fig, ax
