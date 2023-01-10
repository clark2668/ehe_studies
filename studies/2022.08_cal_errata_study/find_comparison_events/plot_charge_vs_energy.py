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

he_file = '/disk20/users/brian/IceCube/juliet/numu_high_energy_merged_999files.hdf5'
he_file = '/disk20/users/brian/IceCube/juliet/nue_high_energy_merged_998files.hdf5'

with tables.open_file(he_file) as f:

    charge = f.get_node('/Homogenized_QTot')
    charge = charge.col('value')    

    primary = f.get_node('/I3JulietPrimaryParticle')
    energies = primary.col('energy')
    
    weight_dict, prop_matrix, evts_per_file = weighting.get_juliet_weightdict_and_propmatrix(f)

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
    
    def GeV_to_PE(e_in_gev, conversion=0.075):
        # rough conversion from HESE is 100 PE = 1 TeV
        # so 0.1 PE per GeV
        return conversion * e_in_gev

    
    e_plot_bins = np.linspace(2, 8, 25)
    charge_plot_bins = GeV_to_PE(np.power(10., e_plot_bins))
    charge_cut_bins = GeV_to_PE(np.power(10., e_plot_bins), 0.025)

    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    bins = [np.linspace(4,10,60), e_plot_bins]
    vals, x_edges, y_edges, im = ax2.hist2d(
        np.log10(energies), np.log10(charge), #weights=weights_from_calc,
        bins=bins, cmap=plt.cm.plasma, norm=colors.LogNorm()
    )
    ax2.plot(e_plot_bins, np.log10(charge_plot_bins), 
             '--', color='red', label='0.075 PE / GeV')
    ax2.plot(e_plot_bins, np.log10(charge_cut_bins), 
             '--', color='blue', label='0.025 PE / GeV')

    ax2.legend()
    sim_cbar = plt.colorbar(im, ax=ax2, label='unweighted Events')
    # lims = 1E-3, 1E-1
    # im.set_clim(lims)
    ax2.set_xlabel('Neutrino Energy log10(GeV)')
    ax2.set_ylabel('log10(HQTot)')
    fig2.tight_layout()
    fig2.savefig('charge_vs_energy.png')
