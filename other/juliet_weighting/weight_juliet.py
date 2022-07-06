from email import utils
import h5py
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import pandas as pd
import tables
import weights as weighter

n_gen = 10000 * 150
livetime = 86400 * 365 # 1 year

he_file = 'juliet_weight_test.hd5'
he_file = 'juliet_weight_test_vhe.hd5'

with tables.open_file(he_file) as f:

    primary = f.get_node('/I3JulietPrimaryParticle')
    energies = primary.col('energy')

    weight_dict = f.get_node('/JulietWeightDict')
    prop_matrix = f.get_node('/PropagationMatrixNuTau').col('item').reshape(-1,140)
    prop_matrix += f.get_node('/PropagationMatrixNuMu').col('item').reshape(-1,140)
    prop_matrix += f.get_node('/PropagationMatrixNuE').col('item').reshape(-1,140)

    # print(weight_dict.col('InjectionSurfaceR'))
    # print(weight_dict.colnames)


    # get weights as they are pre-calculated
    flux_key = 'JulietPropagationMatrixNeutrinoFlux24'
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
    evs1,_ ,__ = ax.hist(np.log10(energies), bins=bins, weights=weights_from_dict, histtype='step')
    evs2,_ ,__ = ax.hist(np.log10(energies), bins=bins, weights=weights_from_calc, histtype='step')
    print(evs1/evs2)
    ax.set_yscale('log')
    plt.tight_layout()
    fig.savefig('raw_energy.png')

