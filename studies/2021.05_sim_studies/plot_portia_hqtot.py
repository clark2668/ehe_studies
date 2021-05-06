import h5py
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors


parser = argparse.ArgumentParser()
parser.add_argument("-d", type=str,
	dest="dataset",
	help="what data set are we plotting (really just somethingn we'll ammend to the title")

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
	try:
		hqtot = np.concatenate((hqtot, file_in['Homogenized_QTot']['value']))
		portia = np.concatenate((portia,file_in['EHEPortiaEventSummarySRT']['bestNPEbtw']))
	except:
		print('  Skipping {}'.format(file))
	file_in.close()

# add "1" to keep the log to a strictly positive number
portia = np.asarray(portia)+1
hqtot = np.asarray(hqtot)+1

portia = np.log10(portia)
hqtot = np.log10(hqtot)

# choose start and end point for 
# plot_min = 0.9 * min(np.min(hqtot), np.min(portia))
# plot_max = 1.1 * max(np.max(hqtot), np.max(portia))
plot_min = 2
plot_max = 7

# fig = plt.figure(figsize=(5,5))
# ax = fig.add_subplot(111)
# ax.hist(hqtot, bins=50, histtype='step', label='hqtot')
# ax.hist(portia,bins=50, histtype='step', label='portia')
# ax.set_yscale('log')
# ax.set_xlabel('Charge')
# ax.set_ylabel('Number of Events')
# ax.legend()
# fig.savefig('hist_charge_{}events.png'.format(len(hqtot)))

sizer=12

fig2 = plt.figure(figsize=(9,7))
ax2 = fig2.add_subplot(111)
my_map = plt.cm.plasma
counts, xedges, yedges, im = ax2.hist2d(portia, hqtot, 
		bins=100,
		range = [[plot_min,plot_max],[plot_min,plot_max]],
		cmap=my_map,
		norm=colors.LogNorm(),
		cmin=1
		)
cbar = plt.colorbar(im, ax=ax2)
ax2.set_title('Set {}'.format(args.dataset))
ax2.plot([plot_min,plot_max],[plot_min,plot_max], 'k:')
cbar.set_label('Number of Events', fontsize=sizer)
ax2.set_xlabel(r'log$_{10}$(Portia NPE)', fontsize=sizer)
ax2.set_ylabel(r'log$_{10}$(Homogenized $Q_{tot})$', fontsize=sizer)
ax2.tick_params(labelsize=sizer)
# plt.tight_layout()
fig2.savefig('portianpe_vs_homogenizeqtot_{}.png'.format(args.dataset), 
	edgecolor='none', bbox_inches='tight')
