from re import A
import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import pandas as pd
import tables
from eheanalysis import weighting

n_gen = 10000 * 150
livetime = 86400 * 365 # 1 year

he_file = '/disk20/users/brian/IceCube/juliet/mu_high_energy_merged_1k.hdf5'

with tables.open_file(he_file) as f:

    charge = f.get_node('/Homogenized_QTot')
    charge = charge.col('value')    

    primary = f.get_node('/I3JulietPrimaryParticle')
    energies = primary.col('energy')
    
    weight_dict, prop_matrix, evts_per_file = weighting.get_juliet_weightdict_and_propmatrix(f)

    # get weights as they are pre-calculated by the I3JulietWeigher
    flux_key = 'JulietPropagationMatrixNeutrinoFlux24'
    weights_from_dict = weighting.calc_juliet_weight_from_weight_dict(
        weight_dict = weight_dict,
        flux_key = flux_key,
        n_gen = n_gen,
        livetime = livetime
    )

    # Now, get weights by calculating them ourselves
    # flux 24 is J(E)*E^2 = 10^-8 [GeV/cm2/s/sr]
    def sample_flux(e_in_gev):
        return 1E-8 * (e_in_gev ** -2.)
    
    # convolve over the propagation matrices ourselves
    weights_from_calc = weighting.calc_juliet_weight_from_weight_dict_and_prop_matrix(
        weight_dict = weight_dict,
        prop_matrix = prop_matrix,
        flux = sample_flux,
        n_gen = n_gen,
        livetime = livetime
    )

    fig = plt.figure()
    ax = fig.add_subplot(111)
    bins = np.linspace(4, 12, 60)
    evs1,_ ,__ = ax.hist(np.log10(energies), bins=bins, weights=weights_from_dict, histtype='step', label='Juliet Precalc')
    evs2,_ ,__ = ax.hist(np.log10(energies), bins=bins, weights=weights_from_calc, histtype='step', label='Max/Brian Calculation', linestyle='--')
    ax.set_yscale('log')
    ax.set_xlabel('Primary Energy log10(GeV)')
    ax.set_ylabel('Weighted Counts/Year')
    ax.legend()
    plt.tight_layout()
    fig.savefig('weighted_energy.png')
    del fig, ax

    # fig2 = plt.figure()
    # ax2 = fig2.add_subplot(111)
    # bins = [np.linspace(4,10,60), np.linspace(2, 7, 21)]
    # vals, x_edges, y_edges, im = ax2.hist2d(
    #     np.log10(energies), np.log10(charge), weights=weights_from_calc,
    #     bins=bins, cmap=plt.cm.plasma, norm=colors.LogNorm()
    # )
    # sim_cbar = plt.colorbar(im, ax=ax2, label='Weighted Events')
    # lims = 1E-3, 1E-1
    # im.set_clim(lims)
    # ax2.set_xlabel('Muon Energy at Detector log10(GeV)')
    # ax2.set_ylabel('log10(HQTot)')
    # fig2.savefig('charge_vs_energy.png')
