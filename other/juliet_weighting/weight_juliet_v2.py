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
# he_file = '/home/mmeier/data/simulations/table_based_sim/juliet/mu/high_energy/l2_prime/mu_high_energy_l2_prime_1k_merged.hd5'
he_file = "/home/mmeier/data/simulations/table_based_sim/juliet/mu/high_energy/l2_prime/1/Level2_prime_00000001.hd5"

with tables.open_file(he_file) as f:

    charge = f.get_node('/Homogenized_QTot')
    charge = charge.col('value')
    print(charge)
    
    # primary = f.get_node('/I3JulietPrimaryParticle')
    # energies = primary.col('energy')

    # weight_dict = f.get_node('/JulietWeightDict')
    # prop_matrix = f.get_node('/PropagationMatrixNuTau').col('item').reshape(-1,140)
    # prop_matrix += f.get_node('/PropagationMatrixNuMu').col('item').reshape(-1,140)
    # prop_matrix += f.get_node('/PropagationMatrixNuE').col('item').reshape(-1,140)

    # # get weights by calculating them ourselves
    # def sample_flux(e_in_gev):
    #     return 1E-8 * (e_in_gev ** -2.)
    
    # energy_bins, weights = weighter.calc_juliet_flux_weight(
    #     weight_dict = weight_dict,
    #     prop_matrix = prop_matrix,
    #     flux = sample_flux,
    #     n_gen = n_gen,
    #     livetime = livetime
    # )
    # bin_centers = (energy_bins[1:] + energy_bins[:-1]) * 0.5

    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # evs1,_ ,__ = ax.hist(
    #     bin_centers, 
    #     bins=energy_bins, 
    #     weights=weights, histtype='step'
    #     )
    # ax.set_xscale('log')
    # ax.set_yscale('log')
    # ax.set_xlabel('Neutrino Energy log10(GeV)')
    # ax.set_ylabel('Weighted Counts/Year')
    # plt.tight_layout()
    # fig.savefig('events_vs_enu.png')
    # del fig, ax