import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import pandas as pd
import tables
import weights as weighter


n_gen = 10000 * 150
livetime = 86400 * 365 # 1 year

# he_file = 'juliet_weight_test.hd5'
# he_file = 'juliet_weight_test_vhe.hd5'
# he_file = '/home/mmeier/data/simulations/table_based_sim/juliet/mu/high_energy/l2_prime/1/Level2_prime_00000001.hd5'
he_file = '/home/mmeier/data/simulations/table_based_sim/juliet/mu/high_energy/l2_prime/mu_high_energy_l2_prime_1k_merged.hd5'

with tables.open_file(he_file) as f:

    charge = f.get_node('/Homogenized_QTot')
    charge = charge.col('value')    

    primary = f.get_node('/I3JulietPrimaryParticle')
    energies = primary.col('energy')

    weight_dict = f.get_node('/JulietWeightDict')
    prop_matrix = f.get_node('/PropagationMatrixNuTau').col('item').reshape(-1,140)
    prop_matrix += f.get_node('/PropagationMatrixNuMu').col('item').reshape(-1,140)
    prop_matrix += f.get_node('/PropagationMatrixNuE').col('item').reshape(-1,140)

    # print(weight_dict.colnames)

    # get weights as they are pre-calculated
    flux_key = 'JulietPropagationMatrixNeutrinoFlux4'
    weights_from_dict = weighter.calc_juliet_weight_from_weight_dict(
        weight_dict = weight_dict,
        flux_key = flux_key,
        n_gen = n_gen,
        livetime = livetime
    )

    # get weights by calculating them ourselves
    def sample_flux(e_in_gev):
        return 1E-8 * (e_in_gev ** -2.)
    
    weights_from_calc = weighter.calc_juliet_weight_from_weight_dict_and_prop_matrix(
        weight_dict = weight_dict,
        prop_matrix = prop_matrix,
        flux = sample_flux,
        n_gen = n_gen,
        livetime = livetime
    )

    fig = plt.figure()
    ax = fig.add_subplot(111)
    bins = np.linspace(4, 12, 60)
    # evs1,_ ,__ = ax.hist(np.log10(energies), bins=bins, weights=weights_from_dict, histtype='step')
    evs2,_ ,__ = ax.hist(np.log10(energies), bins=bins, weights=weights_from_calc, histtype='step')
    ax.set_yscale('log')
    ax.set_xlabel('Primary Energy log10(GeV)')
    ax.set_ylabel('Weighted Counts/Year')
    plt.tight_layout()
    fig.savefig('weighted_energy.png')
    del fig, ax

    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    bins = [np.linspace(4,10,60), np.linspace(2, 7, 21)]
    vals, x_edges, y_edges, im = ax2.hist2d(
        np.log10(energies), np.log10(charge), weights=weights_from_calc,
        bins=bins, cmap=plt.cm.plasma, norm=colors.LogNorm()
    )
    sim_cbar = plt.colorbar(im, ax=ax2, label='Weighted Events')
    lims = 1E-3, 1E-1
    im.set_clim(lims)
    ax2.set_xlabel('Muon Energy at Detector log10(GeV)')
    ax2.set_ylabel('log10(HQTot)')
    fig2.savefig('charge_vs_energy.png')
