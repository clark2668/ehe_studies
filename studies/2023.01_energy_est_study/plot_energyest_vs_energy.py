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

he_file = '/disk20/users/brian/IceCube/juliet/millipede_study/numu_high_energy_999_merged.hdf5'

with tables.open_file(he_file) as f:

    which_var='energy'
    if which_var is 'charge':
        node_name = '/Homogenized_QTot'
        label_name = 'HQtot log10(NPE)'
    elif which_var is 'energy': 
        node_name = '/EHEMuMillipede_SplineMPEseed_sum_contained'
        label_name = 'Millipede Unfoldd Energy log10(GeV)'
    
    use_weights = True
        
    energy_est_node = f.get_node(node_name)
    energy_est = energy_est_node.col('value')
    
    primary = f.get_node('/I3JulietPrimaryParticle')
    energies = primary.col('energy')
    
    reco_dir = f.get_node('/EHE_SplineMPE')
    czen = np.cos(reco_dir.col('zenith'))
    
    charge_node = f.get_node('/Homogenized_QTot')
    charge = charge_node.col('value')
    ndoms = f.get_node(f'/CVMultiplicity').col(f'n_hit_doms')
    
    q_mask = charge > 1
    ndoms_mask = ndoms > 100
    q_n_mask = np.logical_and(q_mask, ndoms_mask)

    # Now, get weights by calculating them ourselves
    # flux 24 is J(E)*E^2 = 10^-8 [GeV/cm2/s/sr]
    def sample_flux(e_in_gev):
        return 1E-8 * (e_in_gev ** -2.)
    
    # convolve over the propagation matrices ourselves
    if use_weights:    
        
        weight_dict, prop_matrix, evts_per_file = weighting.get_juliet_weightdict_and_propmatrix(f)
        weights_from_calc = weighting.calc_juliet_weight_from_weight_dict_and_prop_matrix(
            weight_dict = weight_dict,
            prop_matrix = prop_matrix,
            flux = sample_flux,
            n_gen = n_gen,
            livetime = livetime
        )
    else:
        weights_from_calc = np.ones_like(energies)
    
    def GeV_to_PE(e_in_gev, conversion=0.075):
        # rough conversion from HESE is 100 PE = 1 TeV
        # so 0.1 PE per GeV
        return conversion * e_in_gev

    e_plot_bins = np.linspace(2, 9, 30)

    ######
    #### Energy Estimator vs True E
    ######

    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    bins = [np.linspace(4,10,60), e_plot_bins]
    my_map = plt.cm.plasma
    my_norm = colors.LogNorm()
    vals, x_edges, y_edges, im = ax2.hist2d(
        np.log10(energies[q_n_mask]), np.log10(energy_est[q_n_mask]), 
        weights=weights_from_calc[q_n_mask],
        bins=bins, cmap=my_map, norm=my_norm
    )

    sim_cbar = plt.colorbar(im, ax=ax2, label='unweighted events')
    if use_weights:
        lims = 1E-3, 1E-1
        im.set_clim(lims)
        sim_cbar.ax.set_ylabel('weighted events', rotation=270)
        
    ax2.set_xlabel('Neutrino Energy log10(GeV)')
    ax2.set_ylabel(label_name)
    fig2.tight_layout()
    fig2.savefig(f'{which_var}_vs_energy.png')
    del ax2, fig2, im

    ######
    #### Energy Estimator vs cos(zen)
    ######

    # first juliet
    energyest_vs_czen_bins = [ np.linspace(-1, 1, 20), e_plot_bins ]
    fig, ax = plt.subplots(1,1, figsize=(7,5))
    a, b, c, im1 = ax.hist2d(
        czen[q_n_mask], np.log10(energy_est[q_n_mask]), 
        weights=weights_from_calc[q_n_mask],
        bins = energyest_vs_czen_bins, cmap=my_map, norm=my_norm
    )
    cbar = plt.colorbar(im1, ax=ax, label='Evts / Year')
    if use_weights:
        im1.set_clim(lims)
    ax.set_xlabel(r"SplineMPE cos$\theta$")
    ax.set_ylabel(label_name)
    fig.tight_layout()
    fig.savefig(f"{which_var}_vs_czen_juliet.png", dpi=300)
    del fig, ax