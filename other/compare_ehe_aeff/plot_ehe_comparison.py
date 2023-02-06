import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import style

ic40 = np.genfromtxt(
    'ic40.csv', delimiter=',', skip_header=0,
    names=['E', 'Aeff']
)

ic79 = np.genfromtxt(
    'ic79.csv', delimiter=',', skip_header=0,
    names=['E', 'Aeff']
)

ic86 = np.genfromtxt(
    'ic86.csv', delimiter=',', skip_header=0,
    names=['E', 'Aeff']
)

style.use('mpl.style')

fig, (ax, axr) = plt.subplots(2, 1, 
    sharex=True, 
    gridspec_kw={'height_ratios': [3, 1]},
    figsize=(5,7)
    )

ax.plot(ic40['E'], ic40['Aeff'], label='IC40')
ax.plot(ic79['E'], ic79['Aeff'], label='IC79')
ax.plot(ic86['E'], ic86['Aeff'], label='IC86')

ax.set_yscale('log')
ax.set_ylim([1E0, 1E5])
ax.set_xlim([5,11])
ax.legend()
ax.set_ylabel(f'Effective Area [$m^2$]')

axr.plot(ic40['E'], ic40['Aeff']/ic86['Aeff'], label='40/86')
axr.plot(ic79['E'], ic79['Aeff']/ic86['Aeff'], label='79/86')

axr.set_xlabel(f'Neutrino Energy [GeV]')
axr.set_ylabel(f'ICXX/IC86')
axr.set_ylim([0,2])
axr.axhline(1, ls='--', color='k')

fig.tight_layout()
fig.savefig('aeff_comparison.png', dpi=300)