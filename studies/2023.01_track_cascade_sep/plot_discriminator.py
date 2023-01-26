import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import style
import tables
import pandas as pd
import copy
import yaml
import numpy as np

from functools import partial
from ehefluxes import fluxes
from eheanalysis import weighting, plotting
gzk_flux = fluxes.EHEFlux("cosmogenic_ahlers2010_1E18")
gzk_partial = partial(gzk_flux, which_species="nue_sum") 
cr_flux = weighting.get_flux_model('GaisserH4a', 'corsika')
atmo_flux = weighting.get_flux_model('H3a_SIBYLL23C', 'nugen')

style.use('/home/brian/IceCube/ehe/max_tools/EHE_analysis/eheanalysis/ehe.mplstyle')

livetime = 98.07 * 24 * 60 * 60
print(livetime)

cfg_file = 'config.yaml'
cfg_file = yaml.safe_load(open(cfg_file))

juliet_species = ["nue", "mu"]
# juliet_species = ["nue"]
juliet_energy_levels = ["high_energy", "very_high_energy"]

which_method='new'
if which_method is 'new':
    qcut = 27500
    charge_var = 'hqtot'
    classifier_var = 'speed'
    classifier_bins = np.linspace(0, 0.5, 150)
    nue_cumulative_sign = -1
if which_method is 'old':
    qcut = 25000
    charge_var = 'portia'
    classifier_var = 'redchisqu'
    classifier_bins = np.linspace(0, 1000, 601)
    nue_cumulative_sign = 1 


ehe_weights ={
    "nue": np.asarray([]),
    "mu": np.asarray([])
}
ehe_charge = {
    "nue": np.asarray([]),
    "mu": np.asarray([])
}
ehe_classifier = {
    "nue": np.asarray([]),
    "mu": np.asarray([])
}

events_per_file = {
    "nue_high_energy": 600,
    "nue_very_high_energy": 80,
    "mu_high_energy": 150,
    "mu_very_high_energy": 20
}

for s in juliet_species:
    for l in juliet_energy_levels:

        print(f"Working on juliet {s} {l}")

        the_f = tables.open_file(cfg_file['juliet'][s][l]['file'])
        weight_dict, prop_matrix, evts_per_file = weighting.get_juliet_weightdict_and_propmatrix(the_f)

        evts_per_file = events_per_file[f"{s}_{l}"] # override to fix L2 issue

        charge = the_f.get_node(f"/{cfg_file['variables'][charge_var]['variable']}").col(f"{cfg_file['variables'][charge_var]['value']}")
        classifier = the_f.get_node(f"/{cfg_file['variables'][classifier_var]['variable']}").col(f"{cfg_file['variables'][classifier_var]['value']}")
        n_gen = cfg_file['juliet'][s][l]['n_files'] * evts_per_file

        weights = weighting.calc_juliet_weight_from_weight_dict_and_prop_matrix(
            weight_dict=weight_dict, prop_matrix=prop_matrix, 
            flux=gzk_partial, n_gen=n_gen, livetime=livetime,
        )

        ehe_weights[s] = np.concatenate((ehe_weights[s], copy.deepcopy(abs(weights))))
        ehe_charge[s] = np.concatenate((ehe_charge[s], copy.deepcopy(charge)))
        ehe_classifier[s] = np.concatenate((ehe_classifier[s], copy.deepcopy(classifier)))

        the_f.close()


fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15,5))
lw = 2
ehe_mask = {
    "nue": ehe_charge['nue']>qcut,
    "mu": ehe_charge['mu']>qcut
}
ax1.hist( 
    ehe_classifier['nue'][ehe_mask['nue']], bins=classifier_bins, 
    weights=ehe_weights["nue"][ehe_mask['nue']],
    histtype='step', label=r'GZK $\nu_{e}$', linewidth=lw)
n, b, p = ax1.hist( 
    ehe_classifier['mu'][ehe_mask['mu']], bins=classifier_bins, 
    weights=ehe_weights["mu"][ehe_mask['mu']],
    histtype='step', label=r'GZK $\mu$', linewidth=lw)
# ax1.hist( cor_speed, bins=speed_bins, weights=cor_weights,
#         histtype='step', label='Corsika (H4a)', linewidth=lw)
ax1.set_yscale('log')
ax1.set_xlabel(cfg_file['variables'][classifier_var]['label'])
ax1.set_ylabel('Events / {:.2f} days'.format(livetime/(60*60*24)))
ax1.set_ylim([1E-7, 1E5])
ax1.legend()

ax2.hist( 
    ehe_classifier['nue'][ehe_mask['nue']], bins=classifier_bins, 
    weights=ehe_weights["nue"][ehe_mask['nue']],
    histtype='step', linewidth=lw, density=True)
ax2.hist( 
    ehe_classifier['mu'][ehe_mask['mu']], bins=classifier_bins, 
    weights=ehe_weights["mu"][ehe_mask['mu']],
    histtype='step', linewidth=lw, density=True)
# ax2.set_yscale('log')
ax2.set_xlabel(cfg_file['variables'][classifier_var]['label'])
ax2.set_ylabel('Density')

ax3.hist( 
    ehe_classifier['nue'][ehe_mask['nue']], bins=classifier_bins, 
    weights=ehe_weights["nue"][ehe_mask['nue']],
    histtype='step', linewidth=lw, density=True, cumulative=nue_cumulative_sign)
ax3.hist( 
    ehe_classifier['mu'][ehe_mask['mu']], bins=classifier_bins, 
    weights=ehe_weights["mu"][ehe_mask['mu']],
    histtype='step', linewidth=lw, density=True, cumulative=-nue_cumulative_sign)
ax3.set_xlabel(cfg_file['variables'][classifier_var]['label'])
ax3.set_ylabel('CDF')
ax3.axvline(0.27, linestyle='--')

expected_rate = np.sum(ehe_weights['mu'][ehe_mask['mu']])
print(f"Expected mu rate: {expected_rate}")

def scaling(ax):
    ax.set_xlim([0, 0.5])
for ax in [ax1, ax2, ax3]:
    scaling(ax)


fig.tight_layout()
fig.savefig(f'classifier.png')