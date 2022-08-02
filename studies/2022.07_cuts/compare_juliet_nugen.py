'''
Script to compare weighted nugen to weighted juliet

Let's compare nue, so that we hav "good" control

'''

import numpy as np
import tables, yaml
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

juliet_weights = np.asarray([])
juliet_edets = np.asarray([])
juliet_charges = np.asarray([])
juliet_czens = np.asarray([])

juliet_species = ["nue", "numu", "nutau", "mu", "tau"]
juliet_energy_levels = ["high_energy"]
juliet_species = ["nue"]


selection = {
    'flavor': 'nue',
    'juliet_index': 0
}

for s in juliet_species:
    for l in juliet_energy_levels:

        print("Working on {}".format(s))
        the_f = tables.open_file(cfg_file['juliet'][s][l]['file'])
        weight_dict, prop_matrix, evts_per_file = weighting.get_juliet_weightdict_and_propmatrix(the_f)
        prop_matrix = prop_matrix[selection['juliet_index']] # select only the nue prop matrix
        n_gen = cfg_file['juliet'][s][l]['n_files'] * evts_per_file
        edets = the_f.get_node("/I3JulietPrimaryParticle").col("energy")
        q = the_f.get_node(f"/{charge_var}").col(f"{charge_val}")
        czen = np.cos(the_f.get_node(f"/{czen_var}").col(f"{czen_val}"))

        weights = weighting.calc_juliet_weight_from_weight_dict_and_prop_matrix(
            weight_dict=weight_dict, prop_matrix=prop_matrix, 
            flux=astro_flux, n_gen=n_gen, livetime=livetime
        )
        weights *= 2 # juliet wants us to input the nu+nubar flux, so we need to double the weights
        
        juliet_weights = np.concatenate((juliet_weights, weights))
        juliet_edets = np.concatenate((juliet_edets, edets))
        juliet_charges = np.concatenate((juliet_charges, q))
        juliet_czens = np.concatenate((juliet_czens, czen))
        
        the_f.close()


# nugen nue
nugen_file = pd.HDFStore(cfg_file['nugen'][selection['flavor']]['file'])
nugen_weighter = weighting.get_weighter( nugen_file,  'nugen',
    cfg_file['nugen'][selection['flavor']]['n_files']
    )
nugen_edet = nugen_weighter.get_column('PolyplopiaPrimary', 'energy')
nugen_weights = nugen_weighter.get_weights(astro_flux) * livetime
nugen_charges = nugen_weighter.get_column(charge_var, charge_val)
nugen_czens = np.cos(nugen_weighter.get_column(czen_var, czen_val))
nugen_file.close()


# energy at detector (unweighted -- just to see counts)
fig, ax =  plt.subplots(1,1)
bins = np.logspace(3, 9, 60)
n_j, b_j, p_j = ax.hist(juliet_edets, bins=bins, histtype='step', label="Juliet")
n_n, b_n, p_n = ax.hist(nugen_edet, bins=bins,  histtype='step', label="NuGen")
ax.set_ylabel("Events / {:.2f} days".format(livetime/60/60/24))
ax.set_xlabel(r'E$_{det}$ / GeV')
ax.set_yscale('log')
ax.set_xscale('log')
ax.legend(loc='upper left')
ax.set_title('{}'.format(selection['flavor']))
fig.tight_layout()
fig.savefig('plots/hist_compare_edet_unweighted.png')
del fig, ax


# mask energies outside their overlap energy region
nugen_mask = (np.log10(nugen_edet) > 5) & (np.log10(nugen_edet) < 8)
juliet_mask = (np.log10(juliet_edets) > 5) & (np.log10(juliet_edets) < 8)


# energy at detector (weighted)
fig, (ax, axr) =  plt.subplots(2,1, sharex=True, gridspec_kw={'height_ratios': [3,1]})
bins = np.logspace(4, 9, 60)
n_j, b_j, p_j = ax.hist(juliet_edets[juliet_mask], bins=bins, weights=juliet_weights[juliet_mask],
        histtype='step', label="Juliet")
n_n, b_n, p_n = ax.hist(nugen_edet[nugen_mask], bins=bins, weights=nugen_weights[nugen_mask],
        histtype='step', label="NuGen")
ax.set_ylabel("Events / {:.2f} days".format(livetime/60/60/24))
ax.set_yscale('log')
ax.legend()
ax.set_title('{}'.format(selection['flavor']))
plt.setp(ax.get_xticklabels(), visible=False)

bin_centers = plotting.get_bin_centers(bins)
axr.plot(bin_centers, n_j/n_n - 1, 'o')
axr.set_xlabel(r'E$_{det}$ / GeV')
axr.set_ylabel('Juliet/NuGen - 1')
axr.set_xscale('log')
axr.set_ylim([-0.6,0.6])
axr.axhline(0, linestyle='--')

fig.tight_layout()
fig.savefig('plots/hist_compare_edet.png')
del fig, ax, axr

# charge
fig, (ax, axr) =  plt.subplots(2,1, sharex=True, gridspec_kw={'height_ratios': [3,1]})
bins = np.logspace(2, 8, 60)
n_j, b_j, p_j = ax.hist(juliet_charges[juliet_mask], bins=bins, weights=juliet_weights[juliet_mask],
        histtype='step', label="Juliet")
n_n, b_n, p_n = ax.hist(nugen_charges[nugen_mask], bins=bins, weights=nugen_weights[nugen_mask],
        histtype='step', label="NuGen")
ax.set_ylabel("Events / {:.2f} days".format(livetime/60/60/24))
ax.set_yscale('log')
ax.legend()
ax.set_title('{}'.format(selection['flavor']))
plt.setp(ax.get_xticklabels(), visible=False)

bin_centers = plotting.get_bin_centers(bins)
axr.plot(bin_centers, n_j/n_n - 1, 'o')
axr.set_xlabel(r'Charge / PE')
axr.set_ylabel('Juliet/NuGen - 1')
axr.set_xscale('log')
axr.set_ylim([-0.6,0.6])
axr.axhline(0, linestyle='--')

fig.tight_layout()
fig.savefig('plots/hist_compare_charge.png')
del fig, ax, axr

# cosine zenith
fig, (ax, axr) =  plt.subplots(2,1, sharex=True, gridspec_kw={'height_ratios': [3,1]})
bins = np.linspace(-1, 1, 20)
n_j, b_j, p_j = ax.hist(juliet_czens[juliet_mask], bins=bins, weights=juliet_weights[juliet_mask],
        histtype='step', label="Juliet")
n_n, b_n, p_n = ax.hist(nugen_czens[nugen_mask], bins=bins, weights=nugen_weights[nugen_mask],
        histtype='step', label="NuGen")
ax.set_ylabel("Events / {:.2f} days".format(livetime/60/60/24))
ax.set_yscale('log')
ax.legend()
ax.set_title('{}'.format(selection['flavor']))
plt.setp(ax.get_xticklabels(), visible=False)

bin_centers = plotting.get_bin_centers(bins)
axr.plot(bin_centers, n_j/n_n - 1, 'o')
axr.set_xlabel(r'cos($\theta$)')
axr.set_ylabel('Juliet/NuGen - 1')
axr.set_ylim([-0.6,0.6])
axr.axhline(0, linestyle='--')

fig.tight_layout()
fig.savefig('plots/hist_compare_czen.png')
del fig, ax, axr
