'''
Calculate the juliet effective area, including cuts, and save to npz file for plotting later.
'''

import numpy as np
import tables, yaml, copy
from eheanalysis import weighting, plotting, analysis_9yr
import matplotlib.pyplot as plt
from matplotlib import style
import pandas as pd
style.use('/home/brian/IceCube/ehe/max_tools/EHE_analysis/eheanalysis/ehe.mplstyle')

cfg_file = '../config.yaml'
cfg_file = yaml.safe_load(open(cfg_file))

which_cx = 'cteq5'
qmin = 1E3

# juliet
juliet_species = ["nue", "numu", "nutau", "mu", "tau"]
juliet_energy_levels = ["high_energy"]
neutrino_species = ["nue", "numu", "nutau"]

energy_bins, energy_bin_centers = weighting.get_juliet_enu_binning()

nu_areas_sum = {
    "nue": None,
    "numu": None,
    "nutau": None
}
species_areas_sum = {
    
}

for iS, s in enumerate(juliet_species):
    for l in juliet_energy_levels:

        print("Working on juliet {}".format(s))
        the_f = tables.open_file(cfg_file['juliet'][s][l][f'file_{which_cx}'])
        weight_dict, prop_matrix, evts_per_file = weighting.get_juliet_weightdict_and_propmatrix(the_f)
        
        npe = the_f.get_node('/EHEPortiaEventSummarySRT').col('bestNPEbtw')
        fitqual = the_f.get_node('/EHEOpheliaSRT_ImpLF').col('fitQuality')
        recozen = the_f.get_node('/EHEOpheliaParticleSRT_ImpLF').col('zenith')
        energies = the_f.get_node('/I3JulietPrimaryParticle').col('energy')
        truezen = the_f.get_node('/I3JulietPrimaryParticle').col('zenith')
        n_gen = cfg_file['juliet'][s][l]['n_files'] * evts_per_file
        
        muon_bundle_pass = analysis_9yr.muon_bundle_cut_pass_9yr(recozen, npe)
        track_quality_pass = analysis_9yr.track_quality_cut_pass_9yr(fitqual, npe)
        total_pass = np.logical_and(muon_bundle_pass, track_quality_pass)
        muon_bundle_pass_3yr = analysis_9yr.muon_bundle_cut_pass_3yr(recozen, npe)
        total_pass = muon_bundle_pass_3yr

        # first, get the weighting for this species specifically
        area = weighting.calc_juliet_effective_area(
            energies = energies, weight_dict=weight_dict, n_gen = n_gen,
            energy_bins = energy_bins,  prop_matrix=prop_matrix,
            selection_mask = total_pass
        )
        species_areas_sum [s] = copy.deepcopy(area)
        del area

        for iN, n in enumerate(neutrino_species):
            print("  Working on neutrino species {}".format(n))
            prop_i = prop_matrix[iN]
            any_nans = np.isnan(np.asarray(prop_i.col('item').reshape(-1, 140))).any()
            print("    Have any nans {}".format(any_nans))
            area = weighting.calc_juliet_effective_area(
                energies = energies, weight_dict=weight_dict, n_gen = n_gen,
                energy_bins = energy_bins,  prop_matrix=prop_i,
                selection_mask = total_pass
            )
            if nu_areas_sum[n] is None:
                nu_areas_sum[n] = copy.deepcopy(area)
            else:
                nu_areas_sum[n] += copy.deepcopy(area)
            del area
       
        do_cut_plots = False
        if do_cut_plots:
            charge_bins = np.logspace(2, 7, 32)
            czen_bins = np.linspace(-1, 1, 21)
            fitqual_bins = np.linspace(0,300, 100)
            
            fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2, figsize=(10,10))
            ax1.hist2d(
                np.cos(recozen), npe, bins=[czen_bins, charge_bins], cmin=1, 
            )
            ax2.hist2d(
                np.cos(recozen[total_pass]), npe[total_pass], bins=[czen_bins, charge_bins], cmin=1,
            )
            for ax in [ax1, ax2]:
                ax.set_xlabel("cos(zen)")
                ax.set_ylabel('charge/pe')
                ax.set_yscale('log')

            ax3.hist2d(
                fitqual, npe, bins=[fitqual_bins, charge_bins], cmin=1, 
            )
            ax4.hist2d(
                fitqual[total_pass], npe[total_pass], bins=[fitqual_bins, charge_bins], cmin=1,
            )
            for ax in [ax3, ax4]:
                ax.set_xlabel("speed")
                ax.set_ylabel("charge/pe")
                ax.set_yscale('log')
            fig.suptitle(f"{s}")
            plt.tight_layout()
            fig.savefig(f'./hist_before_after_cuts.png')
            del fig, ax1, ax2, ax3, ax4

        
        
        the_f.close()

import pickle as pickle
output = open('juliet_aeffs.pkl', 'wb')
pickle.dump(
    {
        'energy_bins': energy_bins, 
        'energy_bin_centers': energy_bin_centers,
        'juliet_nu_aeffs': nu_areas_sum, 
        'juliet_species_aeffs': species_areas_sum
    },
    output
)
