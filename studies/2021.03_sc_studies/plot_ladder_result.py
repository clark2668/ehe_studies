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



plot_style = 'relative'
if plot_style == 'absolute':
	axs.plot(data_magsix['setting'], data_magsix['truth'], 'C7o', label='Truth')
	

	# # mag six / noble nine comparison
	# axs.plot(data_magsix['setting'], data_magsix['hqtot_all'], 
	# 	'C0-', marker='x', label='HQtot, MagSix, FADC+ATWD')
	# axs.plot(data_nobnine['setting'], data_nobnine['hqtot_all'], 
	# 	'C1--', marker='+', label='HQtot, NobleNine, FADC+ATWD')
	
	# axs.plot(data_magsix['setting'], data_magsix['portia_all'], 
	# 	'C2-.', marker='1', label='Portia, MagSix, FADC+ATWD')
	# axs.plot(data_nobnine['setting'], data_nobnine['portia_all'], 
	# 	'C3:', marker='2', label='Portia, NobleNine, FADC+ATWD')

	# nobnine ATWD+FADC, ATWD, FADC Comparison
	axs.set_title('Noble Nine')
	axs.plot(data_nobnine['setting'], data_nobnine['hqtot_all'], 
		'C0-', marker='x', label='HQtot, FADC+ATWD')
	axs.plot(data_nobnine['setting'], data_nobnine['hqtot_fadc'], 
		'C0-.', marker='+', label='HQtot, FADC')
	axs.plot(data_nobnine['setting'], data_nobnine['hqtot_atwd'], 
		'C0--', marker='1', label='HQtot, ATWD')


	axs.plot(data_nobnine['setting'], data_nobnine['portia_all'], 
		'C1-', marker='x', label='Portia, FADC+ATWD')
	axs.plot(data_nobnine['setting'], data_nobnine['portia_fadc'], 
		'C1-.', marker='+', label='Portia, FADC')
	axs.plot(data_nobnine['setting'], data_nobnine['portia_atwd'], 
		'C1--', marker='1', label='Portia, ATWD')

	axs.set_ylabel('Filter Wheel Transmittivity')


if plot_style=='relative':


	# # mag six / noble nine comparison
	# axs.plot(data_magsix['setting'], data_magsix['hqtot_all']/data_magsix['truth'], 
	# 	'C0-', marker='x', label='HQtot, MagSix, FADC+ATWD')
	# axs.plot(data_nobnine['setting'], data_nobnine['hqtot_all']/data_nobnine['truth'], 
	# 	'C1--', marker='+', label='HQtot, NobleNine, FADC+ATWD')
	
	# axs.plot(data_magsix['setting'], data_magsix['portia_all']/data_magsix['truth'], 
	# 	'C2-.', marker='1', label='Portia, MagSix, FADC+ATWD')
	# axs.plot(data_nobnine['setting'], data_nobnine['portia_all']/data_nobnine['truth'], 
	# 	'C3:', marker='2', label='Portia, NobleNine, FADC+ATWD')


	# nobnine ATWD+FADC, ATWD, FADC Comparison
	axs.set_title('Noble Nine')
	axs.plot(data_nobnine['setting'], data_nobnine['hqtot_all']/data_magsix['truth'], 
		'C0-', marker='x', label='HQtot, FADC+ATWD')
	axs.plot(data_nobnine['setting'], data_nobnine['hqtot_fadc']/data_magsix['truth'], 
		'C0-.', marker='+', label='HQtot, FADC')
	axs.plot(data_nobnine['setting'], data_nobnine['hqtot_atwd']/data_magsix['truth'], 
		'C0--', marker='1', label='HQtot, ATWD')


	axs.plot(data_nobnine['setting'], data_nobnine['portia_all']/data_magsix['truth'], 
		'C1-', marker='x', label='Portia, FADC+ATWD')
	axs.plot(data_nobnine['setting'], data_nobnine['portia_fadc']/data_magsix['truth'], 
		'C1-.', marker='+', label='Portia, FADC')
	axs.plot(data_nobnine['setting'], data_nobnine['portia_atwd']/data_magsix['truth'], 
		'C1--', marker='1', label='Portia, ATWD')

	axs.set_ylabel('Filter Wheel Transmittivity Measured/Predicted')

axs.set_ylim([0,1.1])
axs.set_xlabel('Setting')
axs.set_xticks([2,3,4,5,6])
axs.legend()
plt.tight_layout()
fig.savefig('ladder_plots/result_comparison_{}_digi_compare.png'.format(plot_style), dpi=300)
del fig, axs

