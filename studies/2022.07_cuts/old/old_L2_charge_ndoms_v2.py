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
alt_ndoms_var = 'HitMultiplicityValues'
alt_czen_var = 'LineFit_redo'
do_efficiency = False
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
    return 1.44e-18 / 2 * (energy/1e5)**-2.2

# set up datasets
# juliet (EHE neutrinos)
juliet_species = ["nue", "numu", "nutau", "mu", "tau"]
juliet_energy_levels = ["high_energy"]
# juliet_species = ["nue"]
# juliet_species = ["mu"]
# juliet_species = []

corsika_sets = ["20787"]
# corsika_sets = []

nugen_sets = ["nue", "numu"]
# nugen_sets = []

burn_samples = ["IC86-I-pass2", "IC86-II-pass2", "IC86-III-pass2"]

style.use('/home/brian/IceCube/ehe/max_tools/EHE_analysis/eheanalysis/ehe.mplstyle')

livetime = 0
for b in burn_samples:
    livetime += cfg_file['burn_sample'][b]['livetime']
print("Total livetime {}".format(livetime))

charge_bins = np.logspace(4, 7, 16)
charge_bin_centers = plotting.get_bin_centers(charge_bins)

charge_bins_noqcuts = np.logspace(2, 7, 16)
charge_bin_centers_noqcuts = plotting.get_bin_centers(charge_bins_noqcuts)

ndoms_bins = np.linspace(0,4000, 21)
ndoms_bin_centers = plotting.get_bin_centers(ndoms_bins)

czen_bins = np.linspace(-1, 1, 20)
czen_bin_centers = plotting.get_bin_centers(czen_bins)


#############################
# ehe/cosmogenic flux (juliet)
#############################
ehe_weights = np.asarray([])
ehe_charge = np.asarray([])
ehe_ndoms = np.asarray([])
ehe_czen = np.asarray([])
ehe_L2_mask = np.asarray([])

ehe_energy_weights_nocuts = None
ehe_energy_weights_after_L2 = None

for s in juliet_species:
    for l in juliet_energy_levels:

        print(f"Working on juliet {s}")
        the_f = tables.open_file(cfg_file['juliet'][s][l]['file'])
        weight_dict, prop_matrix, evts_per_file = weighting.get_juliet_weightdict_and_propmatrix(the_f)
        charge = the_f.get_node(f'/{charge_var}').col(f'{charge_val}')
        ndoms = the_f.get_node(f'/{ndoms_var}').col(f'{ndoms_val}')
        czen = np.cos(the_f.get_node(f'/{czen_var}').col(f'{czen_val}'))
        n_gen = cfg_file['juliet'][s][l]['n_files'] * evts_per_file

        weights = weighting.calc_juliet_weight_from_weight_dict_and_prop_matrix(
            weight_dict=weight_dict, prop_matrix=prop_matrix, 
            flux=gzk_partial, n_gen=n_gen, livetime=livetime
        )
        ehe_weights = np.concatenate((ehe_weights, weights))
        ehe_charge = np.concatenate((ehe_charge, charge))
        ehe_ndoms = np.concatenate((ehe_ndoms, ndoms))
        ehe_czen = np.concatenate((ehe_czen, czen))

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
    fig.savefig("plots/eff_L2.png")
    del fig, ax1, ax2, ax3


#############################
# muon bundles (corsika)
#############################
cor_charge = np.asarray([])
cor_ndoms = np.asarray([])
cor_weights = np.asarray([])
cor_czen = np.asarray([])
for c in corsika_sets:
    print("Working on corsika {}".format(c))
    cor_file = pd.HDFStore(cfg_file['corsika'][c]['file'])
    cor_weighter = weighting.get_weighter( cor_file, 'corsika',
        cfg_file['corsika'][c]['n_files']
        )
    this_cor_charge = cor_weighter.get_column(charge_var, charge_val)
    q_mask = this_cor_charge > q_cut_for_plot
    try:
        this_cor_ndoms = cor_weighter.get_column(ndoms_var, ndoms_val)[q_mask]
    except:
        this_cor_ndoms = cor_weighter.get_column(alt_ndoms_var, ndoms_val)[q_mask]

    try:
        this_cor_czen = np.cos(cor_weighter.get_column(czen_var, czen_val))[q_mask]
    except:
        this_cor_czen = np.cos(cor_weighter.get_column(alt_czen_var, czen_val))[q_mask]

    this_cor_weights = cor_weighter.get_weights(cr_flux)[q_mask] * livetime
    
    cor_charge = np.concatenate((cor_charge, this_cor_charge[q_mask]))
    cor_czen = np.concatenate((cor_czen, this_cor_czen))
    cor_ndoms = np.concatenate((cor_ndoms, this_cor_ndoms))
    cor_weights = np.concatenate((cor_weights, this_cor_weights))
    cor_file.close()


#############################
# atmospheric and astrophysical neutrinos (nugen)
#############################
nugen_charges = np.asarray([])
nugen_ndoms = np.asarray([])
nugen_czens = np.asarray([])
nugen_atmo_weights = np.asarray([])
nugen_astro_weights = np.asarray([])
for n in nugen_sets:
    print("Working on nugen {}".format(n))
    nugen_file = pd.HDFStore(cfg_file['nugen'][n]['file'])
    nugen_weighter = weighting.get_weighter( nugen_file,  'nugen',
        cfg_file['nugen'][n]['n_files']
        )
    this_nugen_charge = np.asarray(nugen_weighter.get_column(charge_var, charge_val))
    q_mask = this_nugen_charge > q_cut_for_plot
    
    try:
        this_nugen_ndoms = np.asarray(nugen_weighter.get_column(ndoms_var, ndoms_val))
    except:
        this_nugen_ndoms = np.asarray(nugen_weighter.get_column(alt_ndoms_var, ndoms_val))
    
    try:
        this_nugen_czen = np.cos(nugen_weighter.get_column(czen_var, czen_val))
    except:
        this_nugen_czen = np.cos(cor_weighter.get_column(alt_czen_var, czen_val))
    
    nugen_charges = np.concatenate((nugen_charges, this_nugen_charge[q_mask]))
    nugen_ndoms = np.concatenate((nugen_ndoms, this_nugen_ndoms[q_mask]))
    nugen_czens = np.concatenate((nugen_czens, this_nugen_czen[q_mask]))
    
    this_nugen_atmo_weights = nugen_weighter.get_weights(atmo_flux) * livetime
    nugen_atmo_weights = np.concatenate((nugen_atmo_weights, this_nugen_atmo_weights[q_mask]))

    this_nugen_astro_weights = nugen_weighter.get_weights(astro_flux) * livetime
    nugen_astro_weights = np.concatenate((nugen_astro_weights, this_nugen_astro_weights[q_mask]))
    nugen_file.close()


#############################
# data
#############################
data_charges = np.asarray([])
data_ndoms = np.asarray([])
data_czen = np.asarray([])
for b in burn_samples:
    print("Working on burn samples {}".format(b))
    data_file = pd.HDFStore(cfg_file['burn_sample'][b]['file'])
    this_charge = np.asarray(data_file.get(charge_var).get(charge_val))
    q_mask = this_charge > q_cut_for_plot
    try:
        this_ndoms = np.asarray(data_file.get(ndoms_var).get(ndoms_val))
    except:
        this_ndoms = np.asarray(data_file.get(alt_ndoms_var).get(ndoms_val))
    
    try:
        this_czen = np.cos(np.asarray(data_file.get(czen_var).get(czen_val)))
    except:
        this_czen = np.cos(np.asarray(data_file.get(alt_czen_var).get(czen_val)))

    # mask = (this_charge > q_cut) & (this_ndoms > ndom_cut) & (this_czen < -0.25)    
    # runs = np.asarray(data_file.get("I3EventHeader").get("Run"))
    # events = np.asfarray(data_file.get("I3EventHeader").get("Event"))
    # print(this_czen[mask])
    # print(this_charge[mask])
    # print(runs[mask])
    # print(events[mask])
    
    data_charges = np.concatenate((data_charges, this_charge[q_mask]))
    data_ndoms = np.concatenate((data_ndoms, this_ndoms[q_mask]))
    data_czen = np.concatenate((data_czen, this_czen[q_mask]))
    data_file.close()

data_q_binned, data_bins_q = np.histogram(data_charges, bins=charge_bins)
errs_q = np.sqrt(data_q_binned)

data_ndoms_binned, data_bins_ndoms = np.histogram(data_ndoms, bins=ndoms_bins)
errs_ndoms = np.sqrt(data_ndoms_binned)

if do_plots:

    #############################
    # histogram of charge 
    #############################
    # fig, ax = plt.subplots(1,1)
    fig, (ax, axr) = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [3, 1]})

    ehe_q_mask = ehe_charge > q_cut_for_plot

    # plot the MC (sum, and individual components)
    sim_q_sum, b, p = ax.hist(
        x = np.concatenate((cor_charge, nugen_charges, nugen_charges,ehe_charge[ehe_q_mask])),
        weights = np.concatenate((cor_weights, nugen_atmo_weights, nugen_astro_weights, ehe_weights[ehe_q_mask])),
        bins = charge_bins, histtype='step', label='Sum', linewidth=2)
    ax.hist(cor_charge, bins=charge_bins, weights=cor_weights, 
        histtype='step', label=r'Atm $\mu$ (H4a)', linewidth=2)
    ax.hist(nugen_charges, bins=charge_bins, weights=nugen_atmo_weights, 
        histtype='step', label=r'Atm $\nu$ (H3a, Sibyll 2.3c)', linewidth=2)
    ax.hist(nugen_charges, bins=charge_bins, weights=nugen_astro_weights, 
        histtype='step', label=r'Astro $\nu$ (north tracks, $\nu_{e}$ + $\nu_{\mu}$ only)', linewidth=2)
    ax.hist(ehe_charge[ehe_q_mask], bins=charge_bins, weights=ehe_weights[ehe_q_mask],
        histtype='step', label=r'Cosmo $\nu$ (Ahlers GZK)', linewidth=2)      

    # plot the data
    ax.errorbar(charge_bin_centers, data_q_binned, yerr=errs_q, fmt='ko', label='Burn Sample')

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylabel('Events / {:.2f} days'.format(livetime/(60*60*24)))
    ax.set_ylim([1E-4, 1E5])
    ax.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left", 
        mode="expand", borderaxespad=0, ncol=2
    )

    # and plot data/MC agreement
    ratio_q = data_q_binned/sim_q_sum
    ratio_q_errs = ratio_q * (errs_q/data_q_binned)
    axr.errorbar(charge_bin_centers, ratio_q, yerr=ratio_q_errs, fmt='ko', label='Data/Sim')
    axr.set_xlabel('')
    axr.set_ylabel('Data/Sim')
    axr.set_xlabel('Q / PE')
    axr.set_ylim([0.5, 1.5])
    axr.plot([1E4, 1E7], [1, 1], 'k--') 
    fig.tight_layout()
    fig.savefig('plots/hist_charge.png')


    #############################
    # histogram of ndoms
    #############################

    fig2, (ax2, ax2r) = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [3, 1]})

    # plot the MC (sum, and individual components)
    sim_ndoms_sum, b, p = ax2.hist(
        x = np.concatenate((cor_ndoms, nugen_ndoms, nugen_ndoms,ehe_ndoms[ehe_q_mask])),
        weights = np.concatenate((cor_weights, nugen_atmo_weights, nugen_astro_weights, ehe_weights[ehe_q_mask])),
        bins = ndoms_bins, histtype='step', label='Sum', linewidth=2)
    ax2.hist(cor_ndoms, bins=ndoms_bins, weights=cor_weights, 
        histtype='step', label=r'Atm $\mu$ (H4a)', linewidth=2)
    ax2.hist(nugen_ndoms, bins=ndoms_bins, weights=nugen_atmo_weights, 
        histtype='step', label=r'Atm $\nu$ (H3a, Sibyll 2.3c)', linewidth=2)
    ax2.hist(nugen_ndoms, bins=ndoms_bins, weights=nugen_astro_weights, 
        histtype='step', label=r'Astro $\nu$ (north tracks, $\nu_{e}$ + $\nu_{\mu}$ only)', linewidth=2)
    ax2.hist(ehe_ndoms[ehe_q_mask], bins=ndoms_bins, weights=ehe_weights[ehe_q_mask],
            histtype='step', label=r'Cosmo $\nu$ (Ahlers GZK)', linewidth=2)

    # plot the data
    ax2.errorbar(ndoms_bin_centers, data_ndoms_binned, yerr=errs_ndoms, fmt='ko', label='Burn Sample')

    ax2.set_yscale('log')
    ax2.set_ylabel('Events / {:.2f} days'.format(livetime/(60*60*24)))
    ax2.set_ylim([1E-7, 1E5])
    ax2.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left", 
        mode="expand", borderaxespad=0, ncol=2,
        title="Q > {:.2f}".format(q_cut_for_plot)
    )
    # and plot data/MC agreement
    ratio_ndoms = data_ndoms_binned/sim_ndoms_sum
    ration_ndoms_err = ratio_ndoms * (errs_ndoms/data_ndoms_binned)
    ax2r.errorbar(ndoms_bin_centers, ratio_ndoms, yerr=ration_ndoms_err, fmt='ko', label='Data/Sim')
    ax2r.set_xlabel('')
    ax2r.set_xlabel('N DOMs')
    ax2r.set_ylabel('Data/Sim')
    ax2r.set_ylim([0.5, 1.5])
    ax2r.plot([0, 4000], [1, 1], 'k--')

    fig2.tight_layout()
    fig2.savefig('plots/hist_ndoms.png')



    #############################
    # histogram of cos(zen) 
    #############################
    # fig, ax = plt.subplots(1,1)
    fig3, (ax3, ax3r) = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [3, 1]})

    ehe_L2_mask = (ehe_charge > q_cut) & (ehe_ndoms > ndom_cut)
    cor_L2_mask = (cor_charge > q_cut) & (cor_ndoms > ndom_cut)
    nugen_L2_mask = (nugen_charges > q_cut) & (nugen_ndoms > ndom_cut)

    data_L2_mask = (data_charges > q_cut) & (data_ndoms > ndom_cut)
    data_czen_binned, data_bins_czen = np.histogram(data_czen[data_L2_mask], bins=czen_bins)
    errs_czen = np.sqrt(data_czen_binned)


    # plot the MC (sum, and individual components)
    sim_czen_sum, b, p = ax3.hist(
        x = np.concatenate((cor_czen[cor_L2_mask], nugen_czens[nugen_L2_mask], nugen_czens[nugen_L2_mask],ehe_czen[ehe_L2_mask])),
        weights = np.concatenate((cor_weights[cor_L2_mask], nugen_atmo_weights[nugen_L2_mask], nugen_astro_weights[nugen_L2_mask], ehe_weights[ehe_L2_mask])),
        bins = czen_bins, histtype='step', label='Sum', linewidth=2)
    ax3.hist(cor_czen[cor_L2_mask], bins=czen_bins, weights=cor_weights[cor_L2_mask], 
        histtype='step', label=r'Atm $\mu$ (H4a)', linewidth=2)
    ax3.hist(nugen_czens[nugen_L2_mask], bins=czen_bins, weights=nugen_atmo_weights[nugen_L2_mask], 
        histtype='step', label=r'Atm $\nu$ (H3a, Sibyll 2.3c)', linewidth=2)
    ax3.hist(nugen_czens[nugen_L2_mask], bins=czen_bins, weights=nugen_astro_weights[nugen_L2_mask], 
        histtype='step', label=r'Astro $\nu$ (north tracks, $\nu_{e}$ + $\nu_{\mu}$ only)', linewidth=2)
    ax3.hist(ehe_czen[ehe_L2_mask], bins=czen_bins, weights=ehe_weights[ehe_L2_mask],
        histtype='step', label=r'Cosmo $\nu$ (Ahlers GZK)', linewidth=2)

    # plot the data
    ax3.errorbar(czen_bin_centers, data_czen_binned, yerr=errs_czen, fmt='ko', label='Burn Sample')

    ax3.set_yscale('log')
    ax3.set_ylabel('Events / {:.2f} days'.format(livetime/(60*60*24)))
    ax3.set_ylim([1E-5, 1E4])
    ax3.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left", 
        mode="expand", borderaxespad=0, ncol=2,
        title="L3 (Q > {:d}, NDoms > {:d})".format(int(q_cut), int(ndom_cut))
    )

    # and plot data/MC agreement
    ratio_czen = data_czen_binned/sim_czen_sum
    ratio_czen_errs = ratio_czen * (errs_czen/data_czen_binned)
    ax3r.errorbar(czen_bin_centers, ratio_czen, yerr=ratio_czen_errs, fmt='ko', label='Data/Sim')
    ax3r.set_xlabel('')
    ax3r.set_ylabel('Data/Sim')
    ax3r.set_xlabel(r'cos($\theta$)')
    ax3r.set_ylim([0.5, 1.5])
    ax3r.plot([-1, 1], [1, 1], 'k--')
    fig3.tight_layout()
    fig3.savefig('plots/hist_czen.png')