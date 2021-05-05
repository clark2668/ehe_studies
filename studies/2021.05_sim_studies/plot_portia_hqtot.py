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

hqtot = np.asarray([])
portia = np.asarray([])

for file in files:
	print(file)
	file_in = h5py.File(file, "r")
	hqtot = np.concatenate((hqtot, file_in['Homogenized_QTot']['value']))
	portia = np.concatenate((portia,file_in['EHEPortiaEventSummarySRT']['bestNPEbtw']))
	file_in.close()

portia = np.asarray(portia)
hqtot = np.asarray(hqtot)

portia = np.log10(portia)
hqtot = np.log10(hqtot)

# choose start and end point for 
plot_min = 0.9 * min(np.min(hqtot), np.min(portia))
plot_max = 1.1 * max(np.max(hqtot), np.max(portia))
plot_min = 1
plot_max = 7

fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(111)
ax.hist(hqtot, bins=50, histtype='step', label='hqtot')
ax.hist(portia,bins=50, histtype='step', label='portia')
ax.set_yscale('log')
ax.set_xlabel('Charge')
ax.set_ylabel('Number of Events')
ax.legend()
fig.savefig('hist_charge_{}events.png'.format(len(hqtot)))

fig2 = plt.figure(figsize=(9,7))
ax2 = fig2.add_subplot(111)
my_map = plt.cm.plasma
counts, xedges, yedges, im = ax2.hist2d(portia, hqtot, 
		bins=50,
		range = [[plot_min,plot_max],[plot_min,plot_max]],
		cmap=my_map,
		norm=colors.LogNorm(),
		cmin=1
		)
cbar = plt.colorbar(im, ax=ax2)
cbar.set_label('Number of Events')
ax2.plot([plot_min,plot_max],[plot_min,plot_max], 'k:')

ax2.set_xlabel(r'$log_{10}$(Portia NPE)')
ax2.set_ylabel(r'$log_{10}$(Homogenized $Q_{tot})$')

fig2.savefig('portianpe_vs_homogenizeqtot_{}events.png'.format(len(hqtot)), edgecolor='none', bbox_inches='tight')
