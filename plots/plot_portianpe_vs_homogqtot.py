import h5py
import argparse
import numpy as np
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser()
parser.add_argument("-f", type=str, nargs='+',
	dest="input_files",
	help="paths to input files (absolute or relative)")

args = parser.parse_args()
files = args.input_files

npe = []
homog_qtot = []

for file in files:
	file_in = h5py.File(file, "r")
	data = file_in['data']
	portia_npe = data['portia_npe']
	homogenized_qtot = data['homogenized_qtot']
	for event, this_portia_npe in enumerate(portia_npe):
		npe.append(this_portia_npe)
		homog_qtot.append(homogenized_qtot[event])

npe = np.asarray(npe)
homog_qtot = np.asarray(homog_qtot)

plot_min = np.min(npe)
min_homog_qtot = np.min(homog_qtot)
if(min_homog_qtot>plot_min):
	plot_min = min_homog_qtot
plot_min*=0.9

plot_max = np.max(npe)
max_homog_qtot = np.max(homog_qtot)
if(max_homog_qtot>plot_max):
	plot_max = max_homog_qtot
plot_max*=1.1



fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(111)
ax.hist(npe,bins=50, histtype='step')
ax.set_yscale('log')
ax.set_xlabel('Portia NPE')
ax.set_ylabel('Number of Events')
fig.savefig('test.png')

fig2 = plt.figure(figsize=(9,7))
ax2 = fig2.add_subplot(111)
my_map = plt.cm.plasma
counts, xedges, yedges, im = ax2.hist2d(npe, homog_qtot, 
		bins=50,
		range = [[plot_min,plot_max],[plot_min,plot_max]],
		cmap=my_map,
		cmin=1)
cbar = plt.colorbar(im, ax=ax2)
cbar.set_label('Number of Events')
ax2.plot([plot_min,plot_max],[plot_min,plot_max], 'k:')

# ax2.set_ylim([100,15000])
# ax2.set_xlim([100,15000])

ax2.set_xlabel('Portia NPE')
ax2.set_ylabel(r'Homogenized $Q_{tot}$')

fig2.savefig('portianpe_vs_homogenizeqtot.png', edgecolor='none', bbox_inches='tight')