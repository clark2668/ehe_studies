import matplotlib.pyplot as plt
import numpy as np

data = np.genfromtxt('ladder_results.txt', delimiter=',',
	skip_header=1, names=['setting', 'truth', 'inhq', 'inpo', 'exhq', 'expo'])

fig, axs = plt.subplots(1,1,figsize=(5,5))

plot_style = 'absolute'
if plot_style == 'absolute':
	axs.plot(data['setting'], data['truth'], 'C7o', label='Truth')
	
	axs.plot(data['setting'], data['inhq'], 'C0-', marker='x', label='HQtot, Including ATWD')
	axs.plot(data['setting'], data['exhq'], 'C1--', marker='+', label='HQtot, Excluding ATWD')

	axs.plot(data['setting'], data['inpo'], 'C2-.', marker='1', label='Portia, Including ATWD')
	axs.plot(data['setting'], data['expo'], 'C3:', marker='2', label='Portia, Excluding ATWD')

	axs.set_ylabel('Filter Wheel Transmittivity')


if plot_style=='relative':

	axs.plot(data['setting'], data['inhq']/data['truth'], 'C0-', marker='x', label='HQtot, Including ATWD')
	axs.plot(data['setting'], data['exhq']/data['truth'], 'C1--', marker='+', label='HQtot, Excluding ATWD')

	axs.plot(data['setting'], data['inpo']/data['truth'], 'C2-.', marker='1', label='Portia, Including ATWD')
	axs.plot(data['setting'], data['expo']/data['truth'], 'C3:', marker='2', label='Portia, Excluding ATWD')	

	axs.set_ylabel('Filter Wheel Transmittivity Measured/Predicted')


axs.set_xlabel('Setting')
axs.set_xticks([2,3,4,5,6])
axs.legend()
plt.tight_layout()
fig.savefig('ladder_plots/result_comparison_{}.png'.format(plot_style), dpi=300)
del fig, axs

