import numpy as np
import matplotlib.pyplot as plt
from matplotlib import style
import pickle as pickle
from eheanalysis import weighting, plotting, analysis_9yr
style.use('/home/brian/IceCube/ehe/ehe_software/ehe_code/EHE_analysis/eheanalysis/ehe.mplstyle')

# do this manually...
import sys
sys.path.append('/home/brian/IceCube/ehe/ehe_software/ehe_code/EHE_analysis/neutrino_level/steps/event_selection')
import old_analysis_results as oar 

# some helper functions
from scipy.interpolate import splrep, splev
def spline_aeffs(loges, logaeffs):
    interpolator = splrep(loges, logaeffs)
    return interpolator

def evaluate_splines(interpolator, loges):
    eval_splines = splev(loges, interpolator)
    return eval_splines

   
flavors = ['nue', 'numu', 'nutau']

# get the 9 year values
# analysis_9yr_energies = {
#     'nue': analysis_9yr.nue_es_9yr,
#     'numu': analysis_9yr.numu_es_9yr,
#     'nutau': analysis_9yr.nutau_es_9yr
# }
# analysis_9yr_aeffs = {
#     'nue': analysis_9yr.nue_aeff_9yr,
#     'numu': analysis_9yr.numu_aeff_9yr,
#     'nutau': analysis_9yr.nutau_aeff_9yr
# }

analysis_9yr_energies = {
    'nue': oar.nue_aeff_7yr[:,0],
    'numu': oar.numu_aeff_7yr[:,0],
    'nutau': oar.nutau_aeff_7yr[:,0]
}
analysis_9yr_aeffs = {
    'nue': oar.nue_aeff_7yr[:,1],
    'numu': oar.numu_aeff_7yr[:,1],
    'nutau': oar.nutau_aeff_7yr[:,1]
}

# get the newly computed values
data_juliet = pickle.load(open('juliet_aeffs.pkl', 'br'))
energy_bin_centers = data_juliet['energy_bin_centers']
juliet_nu_aeffs = data_juliet['juliet_nu_aeffs']
for flavor in flavors:
    # collapse onto the e_nu axis
    juliet_nu_aeffs[flavor] = juliet_nu_aeffs[flavor].sum(axis=0)

# figure out the minimum and maximum point of support
# must do per flavor (cuz grr)
lo_es = {}
hi_es = {}

plot_energies_9yr = {}
plot_aeffs_9yr = {}
plot_energies_new = {}
plot_aeffs_new = {}

for flavor in flavors:
    # make sure there is effective volume in the 9yr set
    has_val = analysis_9yr_aeffs[flavor] > 0

    # set up the nine year interpolator
    interpolator_9yr = splrep(
        np.log10(analysis_9yr_energies[flavor][has_val]),
        np.log10(analysis_9yr_aeffs[flavor][has_val])
    )
    print(interpolator_9yr)
    print("------")

#     # find min and max populated in the 9 yr
#     min_9yr = np.min(analysis_9yr_energies[flavor][has_val])
#     max_9yr = np.max(analysis_9yr_energies[flavor][has_val])

#     # make sure there is effective volume in the new data set
#     has_val =  juliet_nu_aeffs[flavor]> 0

#     # find min and max populated in the new set
#     min_new = np.min(energy_bin_centers[has_val])
#     max_new = np.max(energy_bin_centers[has_val])

#     min_global = min([min_9yr, min_new])
#     max_global = max([max_9yr, max_new])

#     # first, mask out the new calculation
#     min_mask = energy_bin_centers > min_global
#     max_mask = energy_bin_centers < max_global
#     total_mask = np.logical_and(min_mask, max_mask)
#     this_energies_new = energy_bin_centers[total_mask]
#     this_aeffs_new = juliet_nu_aeffs[flavor][total_mask]
    
#     # stash these as the things to plot
#     plot_energies_new[flavor] = this_energies_new
#     plot_aeffs_new[flavor] = this_aeffs_new

#     # now, interpolate the 9yr results to this binning
#     this_9yr_aeffs = splev(
#         np.log10(this_energies_new),
#         interpolator_9yr
#     )
#     this_9yr_aeffs = np.power(10., this_9yr_aeffs) # undo the log

#     plot_energies_9yr[flavor] = this_energies_new
#     plot_aeffs_9yr[flavor] = this_9yr_aeffs



# fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2,3, 
#                                                        figsize=(15,10), 
#                                                        sharex=True,
#                                                        gridspec_kw={'height_ratios': [3, 2]})

# linestyles = ['-', '--', ':', '-.', '-', '--']


# # nue comparison
# ax1.plot( plot_energies_new['nue'], plot_aeffs_new['nue'], 
#          linewidth=3, label="New Juliet")
# ax1.plot(plot_energies_9yr['nue'], plot_aeffs_9yr['nue'], 
#          linestyle='--', linewidth=3, label="previous paper")
# ax1.set_title(r"Electron-Type Neutrinos")
# ax4.plot(
#     plot_energies_new['nue'],
#     plot_aeffs_new['nue']/plot_aeffs_9yr['nue'],
#     linewidth=3
# )

# # numu comparison
# ax2.plot( plot_energies_new['numu'], plot_aeffs_new['numu'], 
#          linewidth=3, label="New Juliet")
# ax2.plot(plot_energies_9yr['numu'], plot_aeffs_9yr['numu'], 
#          linestyle='--', linewidth=3, label="previous paper")
# ax2.set_title(r"Muon-Type Neutrinos")
# ax5.plot(
#     plot_energies_new['numu'],
#     plot_aeffs_new['numu']/plot_aeffs_9yr['numu'],
#     linewidth=3
# )

# # nutau comparison
# ax3.plot( plot_energies_new['nutau'], plot_aeffs_new['nutau'], 
#          linewidth=3, label="New Juliet")
# ax3.plot(plot_energies_9yr['nutau'], plot_aeffs_9yr['nutau'], 
#          linestyle='--', linewidth=3, label="previous paper")
# ax3.set_title(r"Tau-Type Neutrinos")
# ax6.plot(
#     plot_energies_new['nutau'],
#     plot_aeffs_new['nutau']/plot_aeffs_9yr['nutau'],
#     linewidth=3
# )

# ax1.legend(loc="upper left", fontsize=12)
# for ax in [ax1, ax2, ax3]:
#     ax.set_xscale('log')
#     ax.set_yscale('log')
#     # ax.set_ylabel(r"Effective Detector Size [$m^2$]")
#     ax.set_ylabel(r"Effective Area [$m^2$]")
#     ax.set_ylim([1E-1, 1E5])
#     ax.set_xlim([1E5, 1E11])

# for ax in [ax4, ax5, ax6]:
#     ax.set_ylim([0.5, 1.5])
#     ax.axhline(1., 0, 1, linestyle='--', color='red')
#     ax.set_ylabel("New/9yr")
#     ax.set_xlabel(r"Neutrino Energy [GeV]")

# plt.tight_layout()
# fig.savefig(f'./eff_area_quant_comparison.png', dpi=300)
# del fig, ax1
