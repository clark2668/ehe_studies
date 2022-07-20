import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import style
import tables
import pandas as pd
import copy
import yaml


cfg_file = 'config.yaml'
cfg_file = yaml.safe_load(open(cfg_file))
charge_var = cfg_file['variables']['charge']['variable']
charge_val = cfg_file['variables']['charge']['value']


# set up our flux models
# for GZK, us partial function, so it's purely a function of energy for weighting
# for this flux model, nue_sum  = numu_sum = nutau_sum, so we can just pick one flux
from functools import partial
from ehefluxes import fluxes
from eheanalysis import weighting, plotting
gzk_flux = fluxes.EHEFlux("ahlers_gzk")
gzk_partial = partial(gzk_flux, which_species="nue_sum") 
cr_flux = weighting.get_flux_model('GaisserH4a', 'corsika')
atmo_flux = weighting.get_flux_model('H3a_SIBYLL23C', 'nugen')

def astro_flux(energy):
    return 1.44e-18 / 2 * (energy/1e5)**-2.2

# set up datasets
# juliet (EHE neutrinos)
juliet_species = ["nue", "numu", "nutau", "mu", "tau"]
juliet_energy_levels = ["high_energy"]
# juliet_species = ["mu"]

corsika_sets = ["20787"]

nugen_sets = ["nue", "numu"]

burn_samples = ["IC86-I-pass2", "IC86-II-pass2", "IC86-III-pass2"]


def get_from_juliet(f, var, val):

    # get the charge
    charge = f.get_node(f'/{var}')
    charge = charge.col(f'{val}')

    # weight dict and prop matrix
    weight_dict = f.get_node('/JulietWeightDict')
    prop_matrix = [f.get_node(f'/PropagationMatrix{flav}') for flav in ['NuE', 'NuMu', 'NuTau']]

    evts_per_file = weight_dict.col('EventsPerFile')[0]

    return charge, weight_dict, prop_matrix, evts_per_file

style.use('/home/brian/IceCube/ehe/max_tools/EHE_analysis/eheanalysis/ehe.mplstyle')

livetime = 0
for b in burn_samples:
    livetime += cfg_file['burn_sample'][b]['livetime']
print("Total livetime {}".format(livetime))

charge_bins = np.logspace(2, 8, 51)
charge_bin_centers = (charge_bins[1:] + charge_bins[:-1]) * 0.5


#############################
# ehe/cosmogenic flux (juliet)
#############################
summed_ehe = None
ehe_q_projection = np.asarray([])
for s in juliet_species:
    for l in juliet_energy_levels:

        print(f"Working on juliet {s}")
        the_f = tables.open_file(cfg_file['juliet'][s][l]['file'])
        charge, weight_dict, prop_matrix, evts_per_file = get_from_juliet(the_f, charge_var, charge_val)
        n_gen = cfg_file['juliet'][s][l]['n_files'] * evts_per_file

        h = weighting.make_enu_2d_hist( weight_dict, n_gen, gzk_partial,
            var_values=charge, var_bins=charge_bins,
            prop_matrix=prop_matrix, livetime=livetime
            )
    
        if summed_ehe is None:
            summed_ehe = copy.deepcopy(h)
        else:
            summed_ehe += copy.deepcopy(h)
    
        the_f.close()
ehe_q_projection = summed_ehe.sum(axis=1) # project down onto the charge axis

#############################
# muon bundles (corsika)
#############################
cor_charge = np.asarray([])
cor_weights = np.asarray([])
for c in corsika_sets:
    print("Working on corsika {}".format(c))
    cor_file = pd.HDFStore(cfg_file['corsika'][c]['file'])
    cor_weighter = weighting.get_weighter( cor_file, 'corsika',
        cfg_file['corsika'][c]['n_files']
        )
    cor_charge = np.concatenate((cor_charge, cor_weighter.get_column(charge_var, charge_val)))
    cor_weights = np.concatenate((cor_weights, cor_weighter.get_weights(cr_flux) * livetime))
    cor_file.close()


#############################
# atmospheric and astrophysical neutrinos (nugen)
#############################
nugen_charges = np.asarray([])
nugen_atmo_weights = np.asarray([])
nugen_astro_weights = np.asarray([])
for n in nugen_sets:
    print("Working on nugen {}".format(n))
    nugen_file = pd.HDFStore(cfg_file['nugen'][n]['file'])
    nugen_weighter = weighting.get_weighter( nugen_file,  'nugen',
        cfg_file['nugen'][n]['n_files']
        )
    this_nugen_charge = np.asarray(nugen_weighter.get_column(charge_var, charge_val))
    nugen_charges = np.concatenate((nugen_charges, this_nugen_charge))
    
    this_nugen_atmo_weights = nugen_weighter.get_weights(atmo_flux) * livetime
    nugen_atmo_weights = np.concatenate((nugen_atmo_weights, this_nugen_atmo_weights))

    this_nugen_astro_weights = nugen_weighter.get_weights(astro_flux) * livetime
    nugen_astro_weights = np.concatenate((nugen_astro_weights, this_nugen_astro_weights))
    nugen_file.close()


#############################
# data
#############################
data_charges = np.asarray([])
for b in burn_samples:
    print("Working on burn samples {}".format(b))
    data_file = pd.HDFStore(cfg_file['burn_sample'][b]['file'])
    this_charge = np.asarray(data_file.get(charge_var).get(charge_val))
    data_charges = np.concatenate((data_charges, this_charge))
    data_file.close()

data, data_bins = np.histogram(data_charges, bins=charge_bins)
errs = np.sqrt(data)

#############################
# histogram of charge
#############################
fig2, ax2 = plt.subplots(1,1)

# plot the MC (sum, and individual components)
ax2.hist(
    x = np.concatenate((cor_charge, nugen_charges, nugen_charges,charge_bin_centers)),
    weights = np.concatenate((cor_weights, nugen_atmo_weights, nugen_astro_weights, ehe_q_projection)),
    bins = charge_bins, histtype='step', label='Sum')
ax2.hist(cor_charge, bins=charge_bins, weights=cor_weights, 
    histtype='step', label=r'Atm $\mu$ (H4a)')
ax2.hist(nugen_charges, bins=charge_bins, weights=nugen_atmo_weights, 
    histtype='step', label=r'Atm $\nu$ (H3a, Sibyll 2.3c)')
ax2.hist(nugen_charges, bins=charge_bins, weights=nugen_astro_weights, 
    histtype='step', label=r'Astro $\nu$ (north tracks, $\nu_{e}$ + $\nu_{\mu}$ only)')
ax2.hist(charge_bin_centers, bins=charge_bins, weights=ehe_q_projection,
    histtype='step', label=r'Cosmo $\nu$ (Ahlers GZK)')

# plot the data
ax2.errorbar(charge_bin_centers, data, yerr=errs, fmt='ko', label='Burn Sample')

ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_ylabel('Events / bin/ {:.2f} days'.format(livetime/(60*60*24)))
ax2.set_xlabel('Q / PE')
ax2.set_ylim([1E-4, 1E7])
ax2.legend()
fig2.tight_layout()
fig2.savefig('q_hist.png')