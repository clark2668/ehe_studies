import matplotlib.pyplot as plt
import numpy as np

data_magsix = np.genfromtxt('ladder_results_magsix.csv', delimiter=',',
	skip_header=1, names=['setting', 'truth', 
							'hqtot_all', 'hqtot_atwd', 'hqtot_fadc', 
							'portia_all', 'portia_atwd', 'portia_fadc'])

data_nobnine = np.genfromtxt('ladder_results_nobnine.csv', delimiter=',',
	skip_header=1, names=['setting', 'truth', 
							'hqtot_all', 'hqtot_atwd', 'hqtot_fadc', 
							'portia_all', 'portia_atwd', 'portia_fadc'])


fig, axs = plt.subplots(1,1,figsize=(5,5))

plot_style = 'absolute'
if plot_style == 'absolute':
	axs.plot(data_magsix['setting'], data_magsix['truth'], 'C7o', label='Truth')
	
	
	# magsix plotting
	# axs.plot(data['setting'], data['inhq'], 'C0-', marker='x', label='HQtot, Including ATWD')
	# axs.plot(data['setting'], data['exhq'], 'C1--', marker='+', label='HQtot, Excluding ATWD')
	# axs.plot(data['setting'], data['inpo'], 'C2-.', marker='1', label='Portia, Including ATWD')
	# axs.plot(data['setting'], data['expo'], 'C3:', marker='2', label='Portia, Excluding ATWD')

	axs.plot(data_magsix['setting'], data_magsix['hqtot_all'], 
		'C0-', marker='x', label='HQtot, MagSix, FADC+ATWD')
	axs.plot(data_nobnine['setting'], data_nobnine['hqtot_all'], 
		'C1--', marker='+', label='HQtot, NobleNine, FADC+ATWD')
	
	axs.plot(data_magsix['setting'], data_magsix['portia_all'], 
		'C2-.', marker='1', label='Portia, MagSix, FADC+ATWD')
	axs.plot(data_nobnine['setting'], data_nobnine['portia_all'], 
		'C3:', marker='2', label='Portia, NobleNine, FADC+ATWD')

	# nobnine plotting

	axs.set_ylabel('Filter Wheel Transmittivity')


if plot_style=='relative':

	# magsix plotting
	# axs.plot(data['setting'], data['inhq']/data['truth'], 'C0-', marker='x', label='HQtot, Including ATWD')
	# axs.plot(data['setting'], data['exhq']/data['truth'], 'C1--', marker='+', label='HQtot, Excluding ATWD')
	# axs.plot(data['setting'], data['inpo']/data['truth'], 'C2-.', marker='1', label='Portia, Including ATWD')
	# axs.plot(data['setting'], data['expo']/data['truth'], 'C3:', marker='2', label='Portia, Excluding ATWD')	

	axs.plot(data_magsix['setting'], data_magsix['hqtot_all']/data_magsix['truth'], 
		'C0-', marker='x', label='HQtot, MagSix, FADC+ATWD')
	axs.plot(data_nobnine['setting'], data_nobnine['hqtot_all']/data_nobnine['truth'], 
		'C1--', marker='+', label='HQtot, NobleNine, FADC+ATWD')
	
	axs.plot(data_magsix['setting'], data_magsix['portia_all']/data_magsix['truth'], 
		'C2-.', marker='1', label='Portia, MagSix, FADC+ATWD')
	axs.plot(data_nobnine['setting'], data_nobnine['portia_all']/data_nobnine['truth'], 
		'C3:', marker='2', label='Portia, NobleNine, FADC+ATWD')

	axs.set_ylabel('Filter Wheel Transmittivity Measured/Predicted')

axs.set_ylim([0,1.1])
axs.set_xlabel('Setting')
axs.set_xticks([2,3,4,5,6])
axs.legend()
plt.tight_layout()
fig.savefig('ladder_plots/result_comparison_{}.png'.format(plot_style), dpi=300)
del fig, axs

