import numpy as np
import tables, yaml, copy
from eheanalysis import weighting, plotting
import matplotlib.pyplot as plt
from matplotlib import style
import pandas as pd
style.use('/home/brian/IceCube/ehe/max_tools/EHE_analysis/eheanalysis/ehe.mplstyle')

cfg_file = 'config.yaml'
cfg_file = yaml.safe_load(open(cfg_file))
charge_var = cfg_file['variables']['charge']['variable']
charge_val = cfg_file['variables']['charge']['value']
czen_var = cfg_file['variables']['zenith']['variable']
czen_val = cfg_file['variables']['zenith']['value']


livetime = 60*60*24*365 # 1 year

def astro_flux(energy):
    # flux of mu @ 100 TeV (basically the per-flavor, per-particle flux)
    return 1.44e-18 / 2 * (energy/1e5)**-2.37


# juliet


juliet_species = ["nue", "numu", "nutau", "mu", "tau"]
juliet_energy_levels = ["high_energy"]
juliet_species = ["nue"]

which_cx = 'cteq5'
qmin = 1E3

selection = {
    'flavor': 'nue',
    'juliet_index': 0
}

juliet_czentrue_weights = None
juliet_charge_weights = None

czen_bins = np.linspace(-1.1, 1.1, 20)
charge_bins = np.logspace(2, 7, 32)

for s in juliet_species:
    for l in juliet_energy_levels:

        print("Working on juliet {}".format(s))
        the_f = tables.open_file(cfg_file['juliet'][s][l][f'file_{which_cx}'])
        weight_dict, prop_matrix, evts_per_file = weighting.get_juliet_weightdict_and_propmatrix(the_f)
        prop_matrix = prop_matrix[selection['juliet_index']] # select only the nue prop matrix
        n_gen = cfg_file['juliet'][s][l]['n_files'] * evts_per_file
        edets = the_f.get_node("/I3JulietPrimaryParticle").col("energy")
        charge = the_f.get_node(f"/{charge_var}").col(f"{charge_val}")
        czen_true = np.cos(the_f.get_node("/I3JulietPrimaryParticle").col("zenith"))

        juliet_mask = (charge > qmin)
        weights_czentrue = weighting.make_enu_2d_hist( weight_dict, n_gen, astro_flux,
            var_values=czen_true, var_bins=czen_bins,
            prop_matrix=prop_matrix, livetime=livetime,
            selection_mask = juliet_mask
            )
        weights_charge = weighting.make_enu_2d_hist( weight_dict, n_gen, astro_flux,
            var_values=charge, var_bins=charge_bins,
            prop_matrix=prop_matrix, livetime=livetime,
            selection_mask = juliet_mask
            )
        weights_czentrue *= 2
        weights_charge *= 2
        if juliet_czentrue_weights is None:
            juliet_czentrue_weights = copy.deepcopy(weights_czentrue)
            juliet_charge_weights = copy.deepcopy(weights_charge)
        else:
            juliet_czentrue_weights += copy.deepcopy(weights_czentrue)
            juliet_charge_weights += copy.deepcopy(weights_charge)        
        the_f.close()

import juliet_nugen_helper as jnh
juliet_enu_weights = juliet_czentrue_weights.sum(axis=0) # project onto enu axis
juliet_czentrue_weights = jnh.mask_and_collapse(juliet_czentrue_weights)
juliet_charge_weights = jnh.mask_and_collapse(juliet_charge_weights)

# nugen nue
print(f"Working on nugen {selection['flavor']}")
nugen_file = pd.HDFStore(cfg_file['nugen'][selection['flavor']]['file'])
nugen_weighter = weighting.get_weighter( nugen_file,  'nugen',
    cfg_file['nugen'][selection['flavor']]['n_files']
    )
nugen_enu = nugen_weighter.get_column('PolyplopiaPrimary', 'energy')
nugen_weights = nugen_weighter.get_weights(astro_flux) * livetime
nugen_charges = nugen_weighter.get_column(charge_var, charge_val)
# nugen_czens = np.cos(nugen_weighter.get_column(czen_var, czen_val))
nugen_czenstrue = np.cos(nugen_weighter.get_column('PolyplopiaPrimary', 'zenith'))
nugen_file.close()

nugen_mask = (np.log10(nugen_enu) > 5.5) & (np.log10(nugen_enu) < 7.5) & (nugen_charges > qmin)


####################
# true cosine zenith
####################
fig, ax, axr = jnh.make_1D_juliet_nugen_comparison(
    bins=czen_bins, juliet_weights = juliet_czentrue_weights,
    nugen_data = nugen_czenstrue[nugen_mask], nugen_weights=nugen_weights[nugen_mask]
)
ax.legend(loc="upper left")
ax.set_ylabel("Events / {:.2f} days".format(livetime/60/60/24))
ax.set_title('{}, juliet with {}'.format(selection['flavor'], which_cx))
axr.set_xlabel(r'cos($\theta$)')
axr.set_ylabel('Juliet/NuGen - 1')
axr.set_ylim([-0.25,0.25])
fig.tight_layout()
fig.savefig(f'plots/hist_compare_czentrue_{which_cx}.png', dpi=300)
del fig, ax, axr

####################
# reco charge
####################
fig, ax, axr = jnh.make_1D_juliet_nugen_comparison(
    bins=charge_bins, juliet_weights = juliet_charge_weights,
    nugen_data = nugen_charges[nugen_mask], nugen_weights=nugen_weights[nugen_mask]
)
ax.legend(loc="upper left")
ax.set_yscale('log')
ax.set_title('{}, juliet with {}'.format(selection['flavor'], which_cx))
ax.set_ylabel("Events / {:.2f} days".format(livetime/60/60/24))
ax.set_ylim([1E-2, 1E0])
axr.set_xscale('log')
axr.set_xlabel(r'Charge / PE')
axr.set_ylabel('Juliet/NuGen - 1')
axr.set_ylim([-0.4,0.4])
fig.tight_layout()
fig.savefig(f'plots/hist_compare_charge_{which_cx}.png', dpi=300)
del fig, ax, axr

nugen_mask = (nugen_charges > qmin)

####################
# weighted surface fluxes vs enu
####################
juliet_e_bins, juliet_e_bins_centers = weighting.get_juliet_enu_binning()
fig, ax, axr = jnh.make_1D_juliet_nugen_comparison(
    bins=juliet_e_bins, juliet_weights = juliet_enu_weights,
    nugen_data = nugen_enu[nugen_mask], nugen_weights=nugen_weights[nugen_mask]
)
ax.legend(loc="upper left")
ax.set_yscale('log')
ax.set_title('{}, juliet with {}'.format(selection['flavor'], which_cx))
ax.set_ylabel("Events / {:.2f} days".format(livetime/60/60/24))
ax.set_ylim([1E-3, 1E1])
axr.set_xscale('log')
axr.set_xlabel(r'E$_{\nu}$ / GeV')
axr.set_ylabel('Juliet/NuGen - 1')
axr.set_ylim([-0.35,0.35])
axr.set_xlim([1E5,1E9])
fig.tight_layout()
fig.savefig(f'plots/hist_compare_enu_{which_cx}.png', dpi=300)
del fig, ax, axr


# # energy at detector (unweighted -- just to see counts)
# fig, ax =  plt.subplots(1,1)
# bins = np.logspace(3, 9, 60)
# n_j, b_j, p_j = ax.hist(juliet_edets, bins=bins, histtype='step', label="Juliet")
# n_n, b_n, p_n = ax.hist(nugen_edet, bins=bins,  histtype='step', label="NuGen")
# ax.set_ylabel("Events / {:.2f} days".format(livetime/60/60/24))
# ax.set_xlabel(r'E$_{det}$ / GeV')
# ax.set_yscale('log')
# ax.set_xscale('log')
# ax.legend(loc='upper left')
# ax.set_title('{}, juliet with {}'.format(selection['flavor'], which_cx))
# fig.tight_layout()
# fig.savefig(f'plots/hist_compare_edet_unweighted_{which_cx}.png')
# del fig, ax
