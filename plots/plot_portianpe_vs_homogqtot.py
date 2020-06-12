import h5py
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors


parser = argparse.ArgumentParser()
parser.add_argument("-f", type=str, nargs='+',
	dest="input_files",
	help="paths to input files (absolute or relative)")

args = parser.parse_args()
files = args.input_files

npe = []
homog_qtot = []

for file in files:
	print(file)
	file_in = h5py.File(file, "r")
	data = file_in['data']
	portia_npe = data['portia_npe']
	homogenized_qtot = data['homogenized_qtot']
	for event, this_portia_npe in enumerate(portia_npe):
		npe.append(np.log10(this_portia_npe))
		homog_qtot.append(np.log10(homogenized_qtot[event]))

npe = np.asarray(npe)
homog_qtot = np.asarray(homog_qtot)

# choose start and end point for 
plot_min = 0.9 * min(np.min(homog_qtot), np.min(npe))
plot_max = 1.1 * max(np.max(homog_qtot), np.max(npe))

plot_min = 1
plot_max = 7

fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(111)
ax.hist(npe,bins=50, histtype='step')
ax.set_yscale('log')
ax.set_xlabel('Portia NPE')
ax.set_ylabel('Number of Events')
fig.savefig('portianpe_{}events.png'.format(len(npe)))

fig2 = plt.figure(figsize=(9,7))
ax2 = fig2.add_subplot(111)
my_map = plt.cm.plasma
counts, xedges, yedges, im = ax2.hist2d(npe, homog_qtot, 
		bins=50,
		range = [[plot_min,plot_max],[plot_min,plot_max]],
		cmap=my_map,
		norm=colors.LogNorm(),
		cmin=1)
cbar = plt.colorbar(im, ax=ax2)
cbar.set_label('Number of Events')
ax2.plot([plot_min,plot_max],[plot_min,plot_max], 'k:')

#ax2.set_yscale('log')
#ax2.set_xscale('log')

# ax2.set_ylim([100,15000])
# ax2.set_xlim([100,15000])

ax2.set_xlabel(r'$log_{10}$(Portia NPE)')
ax2.set_ylabel(r'$log_{10}$(Homogenized $Q_{tot})$')

fig2.savefig('portianpe_vs_homogenizeqtot_{}events.png'.format(len(npe)), edgecolor='none', bbox_inches='tight')
