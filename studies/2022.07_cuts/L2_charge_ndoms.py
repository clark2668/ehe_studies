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
ndoms_var =  cfg_file['variables']['ndoms']['variable']
ndoms_val =  cfg_file['variables']['ndoms']['value']
czen_var = cfg_file['variables']['zenith']['variable']
czen_val = cfg_file['variables']['zenith']['value']
speed_var = cfg_file['variables']['speed']['variable']
speed_val = cfg_file['variables']['speed']['value']
do_efficiency = True
do_plots = True

log10_q_cut = np.log10(27500)
q_cut = np.power(10., log10_q_cut)
ndom_cut = 100

q_cut_for_plot = np.power(10., 4)

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
    # flux of mu + mubar @ 100 TeV (basically the per-flavor flux)
    return 1.44e-18 / 2 * (energy/1e5)**-2.37

# set up datasets
# juliet (EHE neutrinos)
juliet_species = ["nue", "numu", "nutau", "mu", "tau"]
juliet_energy_levels = ["high_energy", "very_high_energy"]
# juliet_species = ["nue"]
# juliet_energy_levels = ["high_energy"]


corsika_sets = ["20787"]
# corsika_sets = []

nugen_sets = ["nue", "numu", "nutau"]
nugen_sets = ["nue"]

burn_samples = ["IC86-I-pass2", "IC86-II-pass2", "IC86-III-pass2"]
# burn_samples = ["IC86-II-pass2"]

style.use('/home/brian/IceCube/ehe/ehe_software/ehe_code/EHE_analysis/eheanalysis/ehe.mplstyle')
# plt.rcParams.update({
#     'font.size': '16',
# })

livetime = 0
for b in burn_samples:
    livetime += cfg_file['burn_sample'][b]['livetime']
print("Total livetime {}".format(livetime))

charge_bins = np.logspace(4, 7, 16)
charge_bin_centers = plotting.get_bin_centers(charge_bins)

charge_bins_noqcuts = np.logspace(2, 7, 16)
charge_bin_centers_noqcuts = plotting.get_bin_centers(charge_bins_noqcuts)

ndoms_bins = np.linspace(-10,2000, 21)
ndoms_bin_centers = plotting.get_bin_centers(ndoms_bins)

czen_bins = np.linspace(-1, 1, 21)
czen_bin_centers = plotting.get_bin_centers(czen_bins)

speed_bins = np.linspace(0, 0.5, 50)
speed_bin_centers = plotting.get_bin_centers(speed_bins)


#############################
# ehe/cosmogenic flux (juliet)
#############################
ehe_weights = np.asarray([])
ehe_charge = np.asarray([])
ehe_ndoms = np.asarray([])
ehe_czen = np.asarray([])
ehe_speed = np.asarray([])
ehe_L2_mask = np.asarray([])
ehe_speed = np.asarray([])

ehe_energy_weights_nocuts = None
ehe_energy_weights_after_L2 = None

for s in juliet_species:
    for l in juliet_energy_levels:

        print(f"Working on juliet {s}, {l}")
        the_f = tables.open_file(cfg_file['juliet'][s][l]['file'])
        weight_dict, prop_matrix, evts_per_file = weighting.get_juliet_weightdict_and_propmatrix(the_f)
        charge = the_f.get_node(f'/{charge_var}').col(f'{charge_val}')
        ndoms = the_f.get_node(f'/{ndoms_var}').col(f'{ndoms_val}')
        czen = np.cos(the_f.get_node(f'/{czen_var}').col(f'{czen_val}'))
        speed = the_f.get_node(f'/{speed_var}').col(f'{speed_val}')
        n_gen = cfg_file['juliet'][s][l]['n_files'] * evts_per_file

        weights = weighting.calc_juliet_weight_from_weight_dict_and_prop_matrix(
            weight_dict=weight_dict, prop_matrix=prop_matrix, 
            flux=gzk_partial, n_gen=n_gen, livetime=livetime
        )
        ehe_weights = np.concatenate((ehe_weights, weights))
        ehe_charge = np.concatenate((ehe_charge, charge))
        ehe_ndoms = np.concatenate((ehe_ndoms, ndoms))
        ehe_czen = np.concatenate((ehe_czen, czen))
        ehe_speed = np.concatenate((ehe_speed, speed))

        L2_q_mask = charge > q_cut
        L2_ndoms_mask = ndoms > ndom_cut
        L2_mask = L2_q_mask & L2_ndoms_mask
        
        ehe_L2_mask = np.concatenate((ehe_L2_mask, L2_mask))
                    
        if do_efficiency:
            # if we want to tabulate efficiencies vs energy
            # then we need to project into enu space

            # all events (no cuts)
            the_enu_weights = weighting.calc_juliet_flux_weight(
                weight_dict = weight_dict, prop_matrix = prop_matrix, 
                flux= gzk_partial, n_gen = n_gen, livetime=livetime
            )
            if ehe_energy_weights_nocuts is None:
                ehe_energy_weights_nocuts = copy.deepcopy(the_enu_weights)
            else:
                ehe_energy_weights_nocuts += copy.deepcopy(the_enu_weights)

            # after L2 cuts
            the_enu_weights_L2 = weighting.calc_juliet_flux_weight(
                weight_dict = weight_dict, prop_matrix = prop_matrix, 
                flux= gzk_partial, n_gen = n_gen, livetime=livetime,
                selection_mask = L2_mask
            )
            if ehe_energy_weights_after_L2 is None:
                ehe_energy_weights_after_L2 = copy.deepcopy(the_enu_weights_L2)
            else:
                ehe_energy_weights_after_L2 += copy.deepcopy(the_enu_weights_L2)
        
        the_f.close()

if do_efficiency:

    # eff for the charge and ndom cut separately
    # these are the 1D projections
    ehe_q_binned, ehe_q_bins = np.histogram(ehe_charge, 
        bins=charge_bins_noqcuts, weights = ehe_weights)
    ehe_ndom_binned, ehe_ndom_bins = np.histogram(ehe_ndoms, 
        bins=ndoms_bins, weights = ehe_weights)

    eff_qcut = plotting.calc_juliet_eff_vs_var(
        ehe_q_binned, charge_bin_centers_noqcuts)
    eff_ndomcut = plotting.calc_juliet_eff_vs_var(
        ehe_ndom_binned, ndoms_bin_centers)

    # eff vs energy    
    energy_bins, energy_bin_centers = weighting.get_juliet_enu_binning()

    total_ehe_weight_nocuts = ehe_energy_weights_nocuts.sum()
    total_ehe_weight_L2 = ehe_energy_weights_after_L2.sum()
    total_eff_L2 = total_ehe_weight_L2/total_ehe_weight_nocuts

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15,5))

    # efficiency for charge and ndom separately
    ax1.plot(charge_bin_centers_noqcuts, eff_qcut, linewidth=3)
    ax1.set_xscale('log')
    ax1.set_xlabel("Charge Cut / PE")
    ax1.set_ylabel("Efficiency")
    ax1.set_ylim([0, 1.1])
    ax1.axvline(q_cut, ls='--')
    
    ax2.plot(ndoms_bin_centers, eff_ndomcut, linewidth=3)
    ax2.set_xlabel("N DOMs Cut")
    ax2.set_ylabel("Efficiency")
    ax2.set_ylim([0, 1.1])
    ax2.axvline(ndom_cut, ls='--')

    # eff vs energy
    ax3.plot(energy_bin_centers, ehe_energy_weights_after_L2/ehe_energy_weights_nocuts,
        linewidth=3
    )
    ax3.set_xscale('log')
    ax3.set_xlabel(r"E$_{\nu}$ / GeV")
    ax3.set_ylabel("Efficiency")
    ax3.set_ylim([0, 1.1])
    title_for_summary = "L2" + r"(Q$_{min}$=" + "{:.0f}, ".format(q_cut) \
        + "NDOM$_{min}$=" + "{:d})".format(ndom_cut) \
        + "\nOverall Eff: {:.2f}%".format(total_eff_L2*100.)
    ax3.set_title(title_for_summary)

    fig.tight_layout()
    fig.savefig("plots/eff_L2.png", dpi=300)
    del fig, ax1, ax2, ax3


# #############################
# # muon bundles (corsika)
# #############################
# cor_charge = np.asarray([])
# cor_ndoms = np.asarray([])
# cor_weights = np.asarray([])
# cor_czen = np.asarray([])
# cor_speed = np.asarray([])
# for c in corsika_sets:
#     print("Working on corsika {}".format(c))
#     cor_file = pd.HDFStore(cfg_file['corsika'][c]['file'])
#     cor_weighter = weighting.get_weighter( cor_file, 'corsika',
#         cfg_file['corsika'][c]['n_files']
#         )
#     this_cor_charge = cor_weighter.get_column(charge_var, charge_val)
#     this_cor_ndoms = cor_weighter.get_column(ndoms_var, ndoms_val)
#     this_cor_czen = np.cos(cor_weighter.get_column(czen_var, czen_val))    
#     this_cor_speed = cor_weighter.get_column(speed_var, speed_val)
#     this_cor_weights = cor_weighter.get_weights(cr_flux) * livetime

#     q_mask = this_cor_charge > q_cut_for_plot
#     cor_charge = np.concatenate((cor_charge, this_cor_charge[q_mask]))
#     cor_czen = np.concatenate((cor_czen, this_cor_czen[q_mask]))
#     cor_ndoms = np.concatenate((cor_ndoms, this_cor_ndoms[q_mask]))
#     cor_weights = np.concatenate((cor_weights, this_cor_weights[q_mask]))
#     cor_speed = np.concatenate((cor_speed, this_cor_speed[q_mask]))
#     cor_file.close()


# #############################
# # atmospheric and astrophysical neutrinos (nugen)
# #############################
# nugen_charges = np.asarray([])
# nugen_ndoms = np.asarray([])
# nugen_czens = np.asarray([])
# nugen_speeds = np.asarray([])
# nugen_atmo_weights = np.asarray([])
# nugen_astro_weights = np.asarray([])
# for n in nugen_sets:
#     print("Working on nugen {}".format(n))
#     nugen_file = pd.HDFStore(cfg_file['nugen'][n]['file'])
#     nugen_weighter = weighting.get_weighter( nugen_file,  'nugen',
#         cfg_file['nugen'][n]['n_files']
#         )
#     this_nugen_charge = np.asarray(nugen_weighter.get_column(charge_var, charge_val))
#     q_mask = this_nugen_charge > q_cut_for_plot
    
#     this_nugen_ndoms = np.asarray(nugen_weighter.get_column(ndoms_var, ndoms_val))
#     this_nugen_czen = np.cos(nugen_weighter.get_column(czen_var, czen_val))
#     this_nugen_speed = nugen_weighter.get_column(speed_var, speed_val)
       
#     nugen_speeds = np.concatenate((nugen_speeds, this_nugen_speed[q_mask]))
#     nugen_charges = np.concatenate((nugen_charges, this_nugen_charge[q_mask]))
#     nugen_ndoms = np.concatenate((nugen_ndoms, this_nugen_ndoms[q_mask]))
#     nugen_czens = np.concatenate((nugen_czens, this_nugen_czen[q_mask]))
    
#     this_nugen_atmo_weights = nugen_weighter.get_weights(atmo_flux) * livetime
#     nugen_atmo_weights = np.concatenate((nugen_atmo_weights, this_nugen_atmo_weights[q_mask]))

#     this_nugen_astro_weights = nugen_weighter.get_weights(astro_flux) * livetime
#     nugen_astro_weights = np.concatenate((nugen_astro_weights, this_nugen_astro_weights[q_mask]))
#     nugen_file.close()


# #############################
# # data
# #############################
# data_charges = np.asarray([])
# data_ndoms = np.asarray([])
# data_czen = np.asarray([])
# data_speed = np.asarray([])
# for b in burn_samples:
#     print("Working on burn samples {}".format(b))
#     data_file = pd.HDFStore(cfg_file['burn_sample'][b]['file'])
#     this_charge = np.asarray(data_file.get(charge_var).get(charge_val))
#     q_mask = this_charge > q_cut_for_plot
#     this_ndoms = np.asarray(data_file.get(ndoms_var).get(ndoms_val))
#     this_czen = np.cos(np.asarray(data_file.get(czen_var).get(czen_val)))
#     this_speed = data_file.get(speed_var).get(speed_val)

#     # mask = (this_charge > q_cut) & (this_ndoms > ndom_cut) &  (this_speed < 0.2)
#     # runs = np.asarray(data_file.get("I3EventHeader").get("Run"))
#     # events = np.asfarray(data_file.get("I3EventHeader").get("Event"))
#     # for c, q, r, e, sp in zip(this_czen[mask], this_charge[mask], runs[mask], events[mask], this_speed[mask]):
#     #     print(f"  Run {r}, Event {e}, Charge {q:.2f}, CZen {c:.2f}, Speed {sp:.2f}")
    
#     data_speed = np.concatenate((data_speed, this_speed[q_mask]))
#     data_charges = np.concatenate((data_charges, this_charge[q_mask]))
#     data_ndoms = np.concatenate((data_ndoms, this_ndoms[q_mask]))
#     data_czen = np.concatenate((data_czen, this_czen[q_mask]))
#     data_file.close()

# if do_plots:

#     #############################
#     # histogram of charge 
#     #############################

#     ehe_q_mask = ehe_charge > q_cut_for_plot
#     sim_vars={
#         'cor': cor_charge,
#         'nugen_atmo': nugen_charges,
#         'nugen_astro': nugen_charges,
#         'ehe': ehe_charge[ehe_q_mask]
#     }
#     sim_weights={
#         'cor': cor_weights,
#         'nugen_atmo': nugen_atmo_weights,
#         'nugen_astro': nugen_astro_weights,
#         'ehe': ehe_weights[ehe_q_mask]
#     }
#     sim_labels={
#         'cor': r'Atm $\mu$ (H4a)',
#         'nugen_atmo': r'Atm $\nu$ (H3a, Sibyll 2.3c)',
#         'nugen_astro': r'Astro $\nu$ (north tracks, $\nu_{e}$ + $\nu_{\mu}$ only)',
#         'ehe': r'Cosmo $\nu$ (Ahlers GZK)'
#     }
    
#     fig, ax, axr = plotting.do_1D_data_mc_comparison(
#         var_bins = charge_bins, sim_vars = sim_vars,
#         sim_weights = sim_weights, sim_labels = sim_labels,
#         data=data_charges
#     )
#     ax.set_xscale('log')
#     ax.set_yscale('log')
#     ax.set_ylim([1E-4, 1E5])
#     ax.set_ylabel('Events / {:.2f} days'.format(livetime/(60*60*24)))
#     axr.set_ylabel('Data/Sim')
#     axr.set_xlabel('Q / PE')
#     ax.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left", 
#         mode="expand", borderaxespad=0, ncol=2)
#     fig.tight_layout()
#     fig.savefig('plots/hist_q.png')
#     del fig, ax, axr
    
    
#     #############################
#     # histogram of ndoms
#     #############################

#     sim_vars={
#         'cor': cor_ndoms,
#         'nugen_atmo': nugen_ndoms,
#         'nugen_astro': nugen_ndoms,
#         'ehe': ehe_ndoms[ehe_q_mask]
#     }
#     fig, ax, axr = plotting.do_1D_data_mc_comparison(
#         var_bins = ndoms_bins, sim_vars = sim_vars,
#         sim_weights = sim_weights, sim_labels = sim_labels,
#         data=data_ndoms
#     )
#     ax.set_yscale('log')
#     ax.set_ylim([1E-7, 1E5])
#     ax.set_ylabel('Events / {:.2f} days'.format(livetime/(60*60*24)))
#     axr.set_xlabel('N DOMs')
#     axr.set_ylabel('Data/Sim')
#     ax.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left", 
#         mode="expand", borderaxespad=0, ncol=2,
#         title="Q > {:.2f}".format(q_cut_for_plot)
#     )
#     fig.tight_layout()
#     fig.savefig('plots/hist_ndoms.png')
#     del fig, ax, axr

    
#     # #############################
#     # # histogram of cos(zen) 
#     # #############################
    
#     # cor_L2_mask = (cor_charge > q_cut) & (cor_ndoms > ndom_cut)
#     # nugen_L2_mask = (nugen_charges > q_cut) & (nugen_ndoms > ndom_cut)
#     # ehe_L2_mask = (ehe_charge > q_cut) & (ehe_ndoms > ndom_cut)
#     # data_L2_mask = (data_charges > q_cut) & (data_ndoms > ndom_cut)

#     # sim_vars={
#     #     'cor': cor_czen[cor_L2_mask],
#     #     'nugen_atmo': nugen_czens[nugen_L2_mask],
#     #     'nugen_astro': nugen_czens[nugen_L2_mask],
#     #     'ehe': ehe_czen[ehe_L2_mask]
#     # }
#     # sim_weights={
#     #     'cor': cor_weights[cor_L2_mask],
#     #     'nugen_atmo': nugen_atmo_weights[nugen_L2_mask],
#     #     'nugen_astro': nugen_astro_weights[nugen_L2_mask],
#     #     'ehe': ehe_weights[ehe_L2_mask]
#     # }
#     # fig, ax, axr = plotting.do_1D_data_mc_comparison(
#     #     var_bins = czen_bins, sim_vars = sim_vars,
#     #     sim_weights = sim_weights, sim_labels = sim_labels,
#     #     data=data_czen[data_L2_mask]
#     # )
#     # ax.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left", 
#     #     mode="expand", borderaxespad=0, ncol=2,
#     #     title="L3 (Q > {:d}, NDoms > {:d})".format(int(q_cut), int(ndom_cut))
#     # )    
#     # ax.set_yscale('log')
#     # ax.set_ylim([1E-5, 1E4])
#     # ax.set_ylabel('Events / {:.2f} days'.format(livetime/(60*60*24)))
#     # axr.set_ylabel('Data/Sim')
#     # axr.set_xlabel(r'cos($\theta$)')
#     # fig.tight_layout()
#     # fig.savefig('plots/hist_czen.png')
#     # del fig, ax, axr


#     #############################
#     # histogram of linespeed
#     #############################
    
#     sim_vars={
#         'cor': cor_speed[cor_L2_mask],
#         'nugen_atmo': nugen_speeds[nugen_L2_mask],
#         'nugen_astro': nugen_speeds[nugen_L2_mask],
#         'ehe': ehe_speed[ehe_L2_mask]
#     }
#     sim_weights={
#         'cor': cor_weights[cor_L2_mask],
#         'nugen_atmo': nugen_atmo_weights[nugen_L2_mask],
#         'nugen_astro': nugen_astro_weights[nugen_L2_mask],
#         'ehe': ehe_weights[ehe_L2_mask]
#     }
#     fig, ax, axr = plotting.do_1D_data_mc_comparison(
#         var_bins = speed_bins, sim_vars = sim_vars,
#         sim_weights = sim_weights, sim_labels = sim_labels,
#         data=data_speed[data_L2_mask]
#     )
#     ax.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left", 
#         mode="expand", borderaxespad=0, ncol=2,
#         title="L3 (Q > {:d}, NDoms > {:d})".format(int(q_cut), int(ndom_cut))
#     )    
#     ax.set_yscale('log')
#     ax.set_ylim([1E-6, 1E5])
#     ax.set_ylabel('Events / {:.2f} days'.format(livetime/(60*60*24)))
#     axr.set_ylabel('Data/Sim')
#     axr.set_xlabel(r'LineFit Speed')
#     fig.tight_layout()
#     fig.savefig('plots/hist_speed.png')
#     del fig, ax, axr

    
#     # finer binning for 2D plots
#     charge_bins = np.logspace(4, 7, 32)
#     charge_bin_centers = plotting.get_bin_centers(charge_bins)
    
#     czen_bins = np.linspace(-1, 1, 41)
#     czen_bin_centers = plotting.get_bin_centers(czen_bins)

#     #############################
#     # 2D: charge vs zenith
#     #############################
    
#     # sim_y_vals = {
#     #     'cor': cor_charge[cor_L2_mask],
#     #     'nugen_atmo': nugen_charges[nugen_L2_mask],
#     #     'nugen_astro': nugen_charges[nugen_L2_mask],
#     #     'ehe': ehe_charge[ehe_L2_mask]
#     # }
#     # sim_x_vals = {
#     #     'cor': cor_czen[cor_L2_mask],
#     #     'nugen_atmo': nugen_czens[nugen_L2_mask],
#     #     'nugen_astro': nugen_czens[nugen_L2_mask],
#     #     'ehe': ehe_czen[ehe_L2_mask]
#     # }
#     # set_labels={
#     #     'cor': r'Atm $\mu$ (H4a)',
#     #     'nugen_atmo': r'Atm $\nu$ (H3a, Sibyll 2.3c)',
#     #     'nugen_astro': r'Astro $\nu$ (north tracks, $\nu_{e}$ + $\nu_{\mu}$ only)',
#     #     'ehe': r'Cosmo $\nu$ (Ahlers GZK)',
#     #     'data': 'Burn Sample ({:.2f} days)'.format(livetime/(60*60*24))
#     # }
#     # data_y_vals = data_charges[data_L2_mask]
#     # data_x_vals = data_czen[data_L2_mask]
#     # my_map = plt.cm.plasma
#     # plotting_options = {
#     #     'xlabel': r'cos($\theta$)',
#     #     'ylabel': 'Q / PE',
#     #     'zlabel': 'Number of Events',
#     #     'cmap': my_map,
#     #     'norm': colors.LogNorm(),
#     #     'zlims':  (1E-5, 1E2)
#     # }
        
#     # fig, plotting_products = plotting.do_2D_data_mc_comparison(
#     #     bins_x=czen_bins, bins_y = charge_bins,
#     #     sim_x_vals=sim_x_vals, sim_y_vals = sim_y_vals,
#     #     sim_weights=sim_weights, set_labels=set_labels,
#     #     data_x_vals=data_x_vals, data_y_vals=data_y_vals,
#     #     plotting_opts = plotting_options
#     # )
#     # for a in plotting_products['axes']:
#     #     plotting_products['axes'][a].set_yscale('log')
#     #     plotting_products['axes'][a].set_xlabel(plotting_options['xlabel'])
#     #     plotting_products['axes'][a].set_ylabel(plotting_options['ylabel'])
#     # for a in plotting_products['cbars']:
#     #     plotting_products['cbars'][a].set_label(plotting_options['zlabel'])
#     # plotting_products['cbars']['ratio'].set_label("Ratio")
        
#     # fig.tight_layout()
#     # fig.savefig('plots/hist2d_q_czen_datamc.png', dpi=300)
#     # del fig, plotting_products    


#     #############################
#     # 2D: charge vs linefit speed
#     #############################
    
#     # sim_y_vals = {
#     #     'cor': cor_charge[cor_L2_mask],
#     #     'nugen_atmo': nugen_charges[nugen_L2_mask],
#     #     'nugen_astro': nugen_charges[nugen_L2_mask],
#     #     'ehe': ehe_charge[ehe_L2_mask]
#     # }
#     # sim_x_vals = {
#     #     'cor': cor_speed[cor_L2_mask],
#     #     'nugen_atmo': nugen_speeds[nugen_L2_mask],
#     #     'nugen_astro': nugen_speeds[nugen_L2_mask],
#     #     'ehe': ehe_speed[ehe_L2_mask]
#     # }
#     # set_labels={
#     #     'cor': r'Atm $\mu$ (H4a)',
#     #     'nugen_atmo': r'Atm $\nu$ (H3a, Sibyll 2.3c)',
#     #     'nugen_astro': r'Astro $\nu$ (north tracks, $\nu_{e}$ + $\nu_{\mu}$ only)',
#     #     'ehe': r'Cosmo $\nu$ (Ahlers GZK)',
#     #     'data': 'Burn Sample ({:.2f} days)'.format(livetime/(60*60*24))
#     # }
#     # data_y_vals = data_charges[data_L2_mask]
#     # data_x_vals = data_speed[data_L2_mask]
#     # my_map = plt.cm.plasma
#     # plotting_options = {
#     #     'xlabel': r'LineFit Speed',
#     #     'ylabel': 'Q / PE',
#     #     'zlabel': 'Number of Events',
#     #     'cmap': my_map,
#     #     'norm': colors.LogNorm(),
#     #     'zlims':  (1E-5, 1E2)
#     # }
        
#     # fig, plotting_products = plotting.do_2D_data_mc_comparison(
#     #     bins_x=speed_bins, bins_y = charge_bins,
#     #     sim_x_vals=sim_x_vals, sim_y_vals = sim_y_vals,
#     #     sim_weights=sim_weights, set_labels=set_labels,
#     #     data_x_vals=data_x_vals, data_y_vals=data_y_vals,
#     #     plotting_opts = plotting_options
#     # )
#     # for a in plotting_products['axes']:
#     #     plotting_products['axes'][a].set_yscale('log')
#     #     plotting_products['axes'][a].set_xlabel(plotting_options['xlabel'])
#     #     plotting_products['axes'][a].set_ylabel(plotting_options['ylabel'])
#     # for a in plotting_products['cbars']:
#     #     plotting_products['cbars'][a].set_label(plotting_options['zlabel'])
#     # plotting_products['cbars']['ratio'].set_label("Ratio")
        
#     # fig.tight_layout()
#     # fig.savefig('plots/hist2d_q_speed_datamc.png', dpi=300)

#     # del fig, plotting_products
