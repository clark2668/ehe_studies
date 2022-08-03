import numpy as np
import tables, yaml, copy
from eheanalysis import weighting, plotting
import matplotlib.pyplot as plt
from matplotlib import style
import pandas as pd

livetime = 60*60*24*365 # 1 year

def astro_flux(energy):
    # flux of mu @ 100 TeV (basically the per-flavor, per-particle flux)
    return 1.44e-18 / 2 * (energy/1e5)**-2.37

czen_bins = np.linspace(-1.1, 1.1, 20)
czen_bins_centers = plotting.get_bin_centers(czen_bins)

##############
# juliet
##############

weights_tot = None
# species_included = ["nue", "numu", "nutau", "mu", "tau"]
species_included = ["nue"]
for s in species_included:

    print(f"Working on juliet {s}")
    the_f = tables.open_file(f"/disk20/users/brian/IceCube/unified_naming/original/{s}_high_energy_l2_prime_1k_merged.hd5")
    weight_dict, prop_matrix, evts_per_file = weighting.get_juliet_weightdict_and_propmatrix(the_f)
    prop_matrix = prop_matrix[0] # select only the nue prop matrix
    n_gen = 1000 * evts_per_file
    juliet_q = the_f.get_node(f"/Homogenized_QTot").col("value")
    juliet_mask = (juliet_q > 1E3)
    juliet_czen = np.cos(the_f.get_node("/I3JulietPrimaryParticle").col("zenith"))

    h_weights = weighting.make_enu_2d_hist( weight_dict, n_gen, astro_flux,
        var_values=juliet_czen, var_bins=czen_bins,
        prop_matrix=prop_matrix, livetime=livetime,
        selection_mask = juliet_mask
        )
    h_weights *= 2 # juliet wants us to input the nu+nubar flux, so we need to double the weights

    if weights_tot is None:
        weights_tot = copy.deepcopy(h_weights)
    else:
        weights_tot += copy.deepcopy(h_weights)
    
    the_f.close()

# set all surface energies outside our region of interest to zero
juliet_e_bins, juliet_e_bins_centers = weighting.get_juliet_enu_binning()
low_mask = np.log10(juliet_e_bins_centers) < 5.5
high_mask = np.log10(juliet_e_bins_centers) > 7.5
juliet_surface_mask = np.logical_or(low_mask, high_mask)
weights_tot[:, juliet_surface_mask] = 0.
juliet_surf_weights = weights_tot.sum(axis=1) # project down onto the czen axis


##############
# nugen nue
##############
print(f"Working on nugen nue ")
import simweights
nugen_file = pd.HDFStore("/disk20/users/brian/IceCube/unified_naming/original/combo_v2_21218.hdf5")
nugen_weighter = simweights.NuGenWeighter(nugen_file, nfiles=11991)
nugen_esurf = nugen_weighter.get_column('PolyplopiaPrimary', 'energy')
nugen_weights = nugen_weighter.get_weights(astro_flux) * livetime
nugen_charges = nugen_weighter.get_column("Homogenized_QTot", "value")
nugen_czens = np.cos(nugen_weighter.get_column('PolyplopiaPrimary', 'zenith'))
nugen_file.close()

nugen_mask = (np.log10(nugen_esurf) > 5.5) & (np.log10(nugen_esurf) < 7.5) & (nugen_charges > 1E3)


##############
# do comparison plot
##############
fig, (ax, axr) =  plt.subplots(2,1, sharex=True, gridspec_kw={'height_ratios': [3,1]})
n_j, b_j, p_j = ax.hist(czen_bins_centers, bins=czen_bins, weights=juliet_surf_weights,
        histtype='step', label="Juliet")
n_n, b_n, p_n = ax.hist(nugen_czens[nugen_mask], bins=czen_bins, weights=nugen_weights[nugen_mask],
        histtype='step', label="NuGen")
ax.set_ylabel("Events / {:.2f} days".format(livetime/60/60/24))
ax.legend(loc="upper left")
ax.set_title("Species Included: {}".format(species_included))
plt.setp(ax.get_xticklabels(), visible=False)

axr.plot(czen_bins_centers, n_j/n_n - 1, 'o')
axr.set_xlabel(r'True cos($\theta$)')
axr.set_ylabel('Juliet/NuGen - 1')
axr.set_ylim([-0.25,0.25])
axr.axhline(0, linestyle='--')

fig.tight_layout()
fig.savefig(f'plots/hist_compare_czen.png')
del fig, ax, axr
