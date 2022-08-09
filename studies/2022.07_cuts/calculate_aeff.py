import numpy as np
import tables, yaml, copy
from eheanalysis import weighting, plotting, cuts
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
juliet_species = ["nue"]

for s in juliet_species:
    for l in juliet_energy_levels:

        print("Working on juliet {}".format(s))
        the_f = tables.open_file(cfg_file['juliet'][s][l][f'file_{which_cx}'])
        weight_dict, prop_matrix, evts_per_file = weighting.get_juliet_weightdict_and_propmatrix(the_f)
        
        npe = the_f.get_node('/EHEPortiaEventSummarySRT').col('bestNPEbtw')
        fitqual = the_f.get_node('/EHEOpheliaSRT_ImpLF').col('fitQuality')
        recozen = the_f.get_node('/EHEOpheliaParticleSRT_ImpLF').col('zenith')
        
        muon_bundle_pass = cuts.original_muon_bundle_cut_pass(recozen, npe)
        track_quality_pass = cuts.original_track_quality_cut_pass(fitqual, npe)
        
        total_pass = np.logical_and(muon_bundle_pass, track_quality_pass)


        # prop_matrix = prop_matrix[selection[flavor_selection]] # select only the relevant prop matrix
        # n_gen = cfg_file['juliet'][s][l]['n_files'] * evts_per_file
        # charge = the_f.get_node(f"/{charge_var}").col(f"{charge_val}")
        # czen_true = np.cos(the_f.get_node("/I3JulietPrimaryParticle").col("zenith"))
        
        charge_bins = np.logspace(2, 7, 32)
        czen_bins = np.linspace(-1, 1, 21)
        fitqual_bins = np.linspace(0,300, 100)

        
        do_cut_plots = False
        if do_cut_plots:
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


