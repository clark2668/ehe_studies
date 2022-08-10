import numpy as np
import tables, yaml, copy
from eheanalysis import weighting, plotting, analysis_9yr
import matplotlib.pyplot as plt
from matplotlib import style
import pandas as pd
style.use('/home/brian/IceCube/ehe/max_tools/EHE_analysis/eheanalysis/ehe.mplstyle')

cfg_file = 'config.yaml'
cfg_file = yaml.safe_load(open(cfg_file))

which_cx = 'cteq5'
qmin = 1E3

# juliet
juliet_species = ["nue", "numu", "nutau", "mu", "tau"]
juliet_energy_levels = ["high_energy"]
# juliet_species = ["numu"]
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

fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2,3, figsize=(15,10))
linestyles = ['-', '--', ':', '-.', '-', '--']
e_mask = energy_bin_centers < 1E9 # block out energies for which we don't have surface fluxes yet

# plot all neutrinos together
for iN, n in enumerate(nu_areas_sum):
    ax1.plot( energy_bin_centers[e_mask], nu_areas_sum[n].sum(axis=0)[e_mask], 
        linestyle=linestyles[iN], label='{}'.format(n), drawstyle='steps-mid', linewidth=2)

# plot all surface fluxes together
for iS, s in enumerate(species_areas_sum):
    ax2.plot( energy_bin_centers[e_mask], species_areas_sum[s].sum(axis=0)[e_mask], 
        linestyle=linestyles[iS], label='{}'.format(s), drawstyle='steps-mid', linewidth=2)

# nue comparison
ax4.plot( energy_bin_centers[e_mask], species_areas_sum['nue'].sum(axis=0)[e_mask], 
         linewidth=2, label="New")
ax4.plot(analysis_9yr.nue_es_9yr, analysis_9yr.nue_aeff_9yr, 
         linestyle='--', linewidth=2, label="2013 paper")
ax4.set_title(r"$\nu_{e}$")

# numu comparison
ax5.plot( energy_bin_centers[e_mask], species_areas_sum['numu'].sum(axis=0)[e_mask], 
         linewidth=2, label="New")
ax5.plot(analysis_9yr.numu_es_9yr, analysis_9yr.numu_aeff_9yr, 
         linestyle='--', linewidth=2, label="2013 paper")
ax5.set_title(r"$\nu_{\mu}$")

# nutau comparison
ax6.plot( energy_bin_centers[e_mask], species_areas_sum['nutau'].sum(axis=0)[e_mask], 
         linewidth=2, label="New")
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
