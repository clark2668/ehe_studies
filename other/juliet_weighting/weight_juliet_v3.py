import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import pandas as pd
import tables
import weights as weighter
import effective_area as ea

plt.rcParams.update({
    'xtick.labelsize': 14, 
    'ytick.labelsize': 14,
    'xtick.major.size': 5, 
    'ytick.major.size': 5,
    'axes.titlesize': 14,
    'axes.labelsize': 14
})

n_gen = 1000 * 150
livetime = 86400 * 365 # 1 year

he_file = '/home/mmeier/data/simulations/table_based_sim/juliet/mu/high_energy/l2_prime/1/Level2_prime_00000001.hd5'
he_file = '/home/mmeier/data/simulations/table_based_sim/juliet/mu/high_energy/l2_prime/mu_high_energy_l2_prime_1k_merged.hd5'

with tables.open_file(he_file) as f:

    # charge = f.get_node('/Homogenized_QTot')
    # charge = charge.col('value')
    charge = f.get_node('/EHEPortiaEventSummarySRT')
    charge = charge.col('bestNPEbtw')
    
    charge_bins = np.logspace(2, 8, 51)
    charge_bin_centers = (charge_bins[1:] + charge_bins[:-1]) * 0.5

    weight_dict = f.get_node('/JulietWeightDict')
    prop_matrix = [f.get_node(f'/PropagationMatrix{flav}') for flav in ['NuE', 'NuMu', 'NuTau']]

    # get weights by calculating them ourselves
    def sample_flux(e_in_gev):
        return 1E-8 * (e_in_gev ** -2.)
    
    gzk_flux = weighter.AhlersGZKFlux()

    enu_bins, enu_bin_centers = ea.get_juliet_enu_binning()

    h = ea.make_enu_2d_hist(
        weight_dict,
        n_gen,
        gzk_flux,
        var_values=charge,
        var_bins=charge_bins,
        prop_matrix=prop_matrix,
        livetime=livetime)

    # make projection along the charge 
    q_projection = h.sum(axis=1)
    q_stepped_path = ea.stepped_path(charge_bins, q_projection)

    # calculate the efficiency
    q_cut = []
    passing_rate = []
    for q in charge_bins:
        mask = charge_bin_centers > q
        passing_rate.append(np.sum(q_projection[mask]))
    q_cut = np.asarray(q_cut)
    passing_rate = np.asarray(passing_rate)

    fig, ax = plt.subplots(1,1)
    X, Y = np.meshgrid(enu_bin_centers, charge_bin_centers)
    im = ax.pcolor(X, Y, h, norm=colors.LogNorm(vmin=1e-8))
    ax.set_xscale('log')
    ax.set_yscale('log')
    cbar = plt.colorbar(im, ax=ax, label='Events / yr')
    ax.set_xlabel(r'$E_{\nu}$ / GeV')
    ax.set_ylabel('Q / PE')
    fig.tight_layout()
    fig.savefig('q_vs_e.png')
    
    fig2, ax2 = plt.subplots(1,1)
    ax2.plot(*q_stepped_path)
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_ylabel(r'Events / yr / bin')
    ax2.set_xlabel('Q / PE')
    fig2.tight_layout()
    fig2.tight_layout()
    fig2.savefig('q_hist.png')
    # ax2.set_ylim([1E-4, 2E-2])

    fig3, ax3 = plt.subplots(1,1)
    ax3.plot(charge_bins, passing_rate/np.sum(q_projection))
    ax3.set_xscale('log')
    ax3.set_ylabel(r'Efficiency')
    ax3.set_xlabel('Q Cut / PE')
    fig3.tight_layout()
    fig3.savefig('q_cut_efficiency.png')

    # prop_matrix = f.get_node('/PropagationMatrixNuTau').col('item').reshape(-1,140)
    # prop_matrix += f.get_node('/PropagationMatrixNuMu').col('item').reshape(-1,140)
    # prop_matrix += f.get_node('/PropagationMatrixNuE').col('item').reshape(-1,140)

    # energy_bins, weights = weighter.calc_juliet_flux_weight(
    #     weight_dict = weight_dict,
    #     prop_matrix = prop_matrix,
    #     flux = sample_flux,
    #     n_gen = n_gen,
    #     livetime = livetime
    # )
    # print(weights)
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
