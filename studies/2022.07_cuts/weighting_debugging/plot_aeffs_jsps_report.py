import numpy as np
import matplotlib.pyplot as plt
from matplotlib import style
from eheanalysis import weighting, plotting, analysis_9yr
style.use('/home/brian/IceCube/ehe/max_tools/EHE_analysis/eheanalysis/ehe.mplstyle')


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


fig, (ax1, ax2, ax3) = plt.subplots(1,3, figsize=(15,5))
linestyles = ['-', '--', ':', '-.', '-', '--']
e_mask = energy_bin_centers < 1E20

# nue comparison
ax1.plot( energy_bin_centers[e_mask], juliet_nu_aeffs['nue'].sum(axis=0)[e_mask], 
        #  linewidth=3, label="fellowship results (this work)")
         linewidth=3, label="New Juliet")
ax1.plot(analysis_9yr.nue_es_9yr, analysis_9yr.nue_aeff_9yr, 
         linestyle='--', linewidth=3, label="previous paper")
ax1.set_title(r"Electron-Type Neutrinos")

# numu comparison
ax2.plot( energy_bin_centers[e_mask], juliet_nu_aeffs['numu'].sum(axis=0)[e_mask], 
         linewidth=3)
ax2.plot(analysis_9yr.numu_es_9yr, analysis_9yr.numu_aeff_9yr, 
         linestyle='--', linewidth=3)
ax2.set_title(r"Muon-Type Neutrinos")

# nutau comparison
ax3.plot( energy_bin_centers[e_mask], juliet_nu_aeffs['nutau'].sum(axis=0)[e_mask], 
         linewidth=3)
ax3.plot(analysis_9yr.nutau_es_9yr, analysis_9yr.nutau_aeff_9yr, 
         linestyle='--', linewidth=3)
ax3.set_title(r"Tau-Type Neutrinos")

plot_nugen = True
if plot_nugen:
    ax1.plot( energy_bin_centers, nugen_nu_aeffs['nue'][0], 
        linewidth=3, linestyle='-.', label='NuGen')
    ax2.plot( energy_bin_centers, nugen_nu_aeffs['numu'][0], 
        linewidth=3, linestyle='-.')
    ax3.plot( energy_bin_centers, nugen_nu_aeffs['nutau'][0], 
        linewidth=3, linestyle='-.')



ax1.legend(loc="upper left", fontsize=12)
for ax in [ax1, ax2, ax3]:
    ax.set_xscale('log')
    ax.set_yscale('log')
    # ax.set_ylabel(r"Effective Detector Size [$m^2$]")
    ax.set_ylabel(r"Effective Area [$m^2$]")
    ax.set_xlabel(r"Neutrino Energy [GeV]")
    ax.set_ylim([1E-1, 1E5])
    ax.set_xlim([1E5, 1E11])

plt.tight_layout()
fig.savefig(f'./eff_area_jsps_report.png', dpi=300)
del fig, ax1
