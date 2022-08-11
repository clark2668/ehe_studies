import numpy as np
import matplotlib.pyplot as plt
from eheanalysis import weighting, plotting, analysis_9yr


import pickle as pickle
data_juliet = pickle.load(open('juliet_aeffs.pkl', 'br'))
energy_bin_centers = data_juliet['energy_bin_centers']
juliet_nu_aeffs = data_juliet['juliet_nu_aeffs']
juliet_species_aeffs = data_juliet['juliet_species_aeffs']

juliet_species = juliet_species_aeffs.keys()
juliet_nu_species = juliet_nu_aeffs.keys()

data_nugen = pickle.load(open('nugen_aeffs.pkl', 'br'))
nugen_nu_aeffs = data_nugen['nugen_nu_aeffs']
nugen_nu_species = nugen_nu_aeffs.keys()


fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2,3, figsize=(15,10))
linestyles = ['-', '--', ':', '-.', '-', '--']
e_mask = energy_bin_centers < 1E9 # block out energies for which we don't have surface fluxes yet

# plot all juliet neutrinos together
for iN, n in enumerate(juliet_nu_species):
    ax1.plot( energy_bin_centers[e_mask], juliet_nu_aeffs[n].sum(axis=0)[e_mask], 
        linestyle=linestyles[iN], label='{}'.format(n), linewidth=2)
ax1.set_title("Juliet")

# plot all juliet surface fluxes together
for iS, s in enumerate(juliet_species):
    ax2.plot( energy_bin_centers[e_mask], juliet_species_aeffs[s].sum(axis=0)[e_mask], 
        linestyle=linestyles[iS], label='{}'.format(s), linewidth=2)
ax1.set_title("Juliet")

for iN, n in enumerate(nugen_nu_species):
    ax3.plot( energy_bin_centers, nugen_nu_aeffs[n][0], 
        linestyle=linestyles[iN], label='{}'.format(n), linewidth=2)
ax3.set_title("NuGen")



# for iN, in in 

# nue comparison
ax4.plot( energy_bin_centers[e_mask], juliet_nu_aeffs['nue'].sum(axis=0)[e_mask], 
         linewidth=2, label="Juliet New")
ax4.plot( energy_bin_centers, nugen_nu_aeffs['nue'][0], 
         linewidth=2, label="NuGen", linestyle='-.')
ax4.plot(analysis_9yr.nue_es_9yr, analysis_9yr.nue_aeff_9yr, 
         linestyle='--', linewidth=2, label="2013 paper")
ax4.set_title(r"$\nu_{e}$")

# numu comparison
ax5.plot( energy_bin_centers[e_mask], juliet_nu_aeffs['numu'].sum(axis=0)[e_mask], 
         linewidth=2, label="Juliet New")
ax5.plot( energy_bin_centers, nugen_nu_aeffs['numu'][0], 
         linewidth=2, label="NuGen", linestyle='-.')
ax5.plot(analysis_9yr.numu_es_9yr, analysis_9yr.numu_aeff_9yr, 
         linestyle='--', linewidth=2, label="2013 paper")
ax5.set_title(r"$\nu_{\mu}$")

# nutau comparison
ax6.plot( energy_bin_centers[e_mask], juliet_nu_aeffs['nutau'].sum(axis=0)[e_mask], 
         linewidth=2, label="Juliet New")
ax6.plot(analysis_9yr.nutau_es_9yr, analysis_9yr.nutau_aeff_9yr, 
         linestyle='--', linewidth=2, label="2013 paper")
ax6.set_title(r"$\nu_{\tau}$")


for ax in [ax1, ax2, ax3, ax4, ax5, ax6]:
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylabel(r"Eff Area / $m^2$")
    ax.set_xlabel(r"E$_{\nu}$ / GeV")
    ax.legend(loc="upper left")
    ax.set_ylim([1E-1, 1E5])
    ax.set_xlim([1E5, 1E11])

plt.tight_layout()
fig.savefig(f'./eff_area.png', dpi=300)
del fig, ax1
