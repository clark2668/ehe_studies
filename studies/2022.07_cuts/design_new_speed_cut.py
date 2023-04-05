'''
Point of this script is to understand how well our proposed
LineFit variable works in Juliet
'''

from sre_parse import fix_flags
import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import style
import tables
import pandas as pd
import copy
import yaml

from eheanalysis import plotting, analysis_9yr, cuts

cfg_file = 'config.yaml'
cfg_file = yaml.safe_load(open(cfg_file))

do_efficiency = True
do_plots = True
version = 'old'
version = 'new'
if version == 'new':
    charge_var = cfg_file['variables']['charge']['variable']
    charge_val = cfg_file['variables']['charge']['value']
    ndoms_var =  cfg_file['variables']['ndoms']['variable']
    ndoms_val =  cfg_file['variables']['ndoms']['value']
    speed_var = cfg_file['variables']['speed']['variable']
    speed_val = cfg_file['variables']['speed']['value']
    speed_bins = np.linspace(0, 2, 251)
    speed_label = "LineFit Speed"
    speed_cut = 0.27
    speed_lims = [0, 0.5]
    log10_q_cut = np.log10(27500)
    nue_cumulative_sign = -1
elif version == 'old':
    charge_var = "EHEPortiaEventSummarySRT"
    charge_val = "bestNPEbtw"
    ndoms_var =  "EHEPortiaEventSummarySRT"
    ndoms_val =  "NCHbtw"
    speed_var = "EHEOpheliaSRT_ImpLF"
    speed_val = "fitQuality"
    speed_label = "Ophelia FitQual"
    speed_bins = np.linspace(0, 1000, 151)
    speed_cut = 100
    speed_lims = [0,500]
    log10_q_cut = np.log10(25000)
    nue_cumulative_sign = 1

charge_bins = np.logspace(4, 7, 31)
charge_bin_centers = plotting.get_bin_centers(charge_bins)
speed_bin_centers = plotting.get_bin_centers(speed_bins)
q_cut = np.power(10., log10_q_cut)
ndom_cut = 100

from functools import partial
from ehefluxes import fluxes
from eheanalysis import weighting
gzk_flux = fluxes.EHEFlux("ahlers_gzk")
gzk_partial = partial(gzk_flux, which_species="nue_sum") 
cr_flux = weighting.get_flux_model('GaisserH4a', 'corsika')
atmo_flux = weighting.get_flux_model('H3a_SIBYLL23C', 'nugen')

burn_samples = ["IC86-I-pass2", "IC86-II-pass2", "IC86-III-pass2"]

style.use('/home/brian/IceCube/ehe/ehe_software/ehe_code/EHE_analysis/eheanalysis/ehe.mplstyle')

livetime = 0
for b in burn_samples:
    livetime += cfg_file['burn_sample'][b]['livetime']
livetime = 203.8 * 24 * 60 * 60 # seconds (203.8 days to seconds)
livetime = 365 * 24 * 60 * 60
print("Total livetime {:1f} ({:.1f} days)".format(livetime, livetime/60/60/24))

#############################
# muon bundles (corsika)
#############################

cor_weights = np.asarray([])
cor_charge = np.asarray([])
cor_speed = np.asarray([])

corsika_sets = ["20787"]
# corsika_sets = []
for c in corsika_sets:
    print("Working on corsika {}".format(c))
    cor_file = pd.HDFStore(cfg_file['corsika'][c]['file'])
    cor_weighter = weighting.get_weighter( cor_file, 'corsika',
        cfg_file['corsika'][c]['n_files']
        )
    this_cor_charge = cor_weighter.get_column(charge_var, charge_val)
    this_cor_ndoms = cor_weighter.get_column(ndoms_var, ndoms_val)
    this_cor_speed = cor_weighter.get_column(speed_var, speed_val)
    
    L2_q_mask = this_cor_charge > q_cut
    L2_ndom_mask = this_cor_ndoms > ndom_cut
    L2_mask = L2_q_mask & L2_ndom_mask # build the L2 cut

    cor_charge = np.concatenate((cor_charge, this_cor_charge[L2_mask]))
    cor_speed = np.concatenate((cor_speed, this_cor_speed[L2_mask]))
    cor_weights = np.concatenate((cor_weights, cor_weighter.get_weights(cr_flux)[L2_mask] * livetime))
    cor_file.close()

#############################
# nugen (for atmospheric neutrino estimates)
#############################

nugen_atmo_weights = np.asarray([])
nugen_charge = np.asarray([])
nugen_speed = np.asarray([])
nugen_energy = np.asarray([])

nugen_sets = ["nue", "numu", "nutau"]
# nugen_sets = ["nue"]
# nugen_sets = []
for n in nugen_sets:
    print("Working on nugen {}".format(n))

    nugen_file = pd.HDFStore(cfg_file['nugen'][n]['file'])
    nugen_weighter = weighting.get_weighter( nugen_file,  'nugen',
        cfg_file['nugen'][n]['n_files']
        )
    this_nugen_charge = np.asarray(nugen_weighter.get_column(charge_var, charge_val))    
    this_nugen_ndoms = np.asarray(nugen_weighter.get_column(ndoms_var, ndoms_val))
    this_nugen_speed = nugen_weighter.get_column(speed_var, speed_val)
    this_nugen_energy = nugen_weighter.get_column('PolyplopiaPrimary', 'energy')
    
    L2_q_mask = this_nugen_charge > q_cut
    L2_ndom_mask = this_nugen_ndoms > ndom_cut
    L2_mask = L2_q_mask & L2_ndom_mask # build the L2 cut
    
    nugen_charge = np.concatenate((nugen_charge, this_nugen_charge[L2_mask]))       
    nugen_speed = np.concatenate((nugen_speed, this_nugen_speed[L2_mask]))
    
    this_nugen_atmo_weights = nugen_weighter.get_weights(atmo_flux) * livetime
    nugen_atmo_weights = np.concatenate((nugen_atmo_weights, this_nugen_atmo_weights[L2_mask]))
    nugen_energy = np.concatenate((nugen_energy, this_nugen_energy[L2_mask]))
    
    nugen_file.close()

summed_nugen_weight = np.sum(nugen_atmo_weights)
print("Expect {} nugen atmo events at L2".format(summed_nugen_weight))

#############################
# ehe/cosmogenic flux (juliet)
#############################

ehe_weights = np.asarray([])
ehe_charge = np.asarray([])
ehe_speed = np.asarray([])
ehe_L2_mask = np.asarray([])

juliet_species = ["nue", "numu", "nutau", "mu", "tau"]
juliet_energy_levels = ["high_energy", "very_high_energy"]
# juliet_species = ["nue"]
# juliet_energy_levels = ["high_energy"]
# juliet_species = []

ehe_energy_weights_nocuts = None
ehe_energy_weights_after_L2 = None
ehe_energy_weights_after_L3 = None

for s in juliet_species:
    for l in juliet_energy_levels:

        if s in juliet_species:

            print(f"Working on juliet {s}, {l}")
            the_f = tables.open_file(cfg_file['juliet'][s][l]['file'])
            weight_dict, prop_matrix, evts_per_file = weighting.get_juliet_weightdict_and_propmatrix(the_f)
            charge = the_f.get_node(f'/{charge_var}').col(f'{charge_val}')
            ndoms = the_f.get_node(f'/{ndoms_var}').col(f'{ndoms_val}')
            speed = the_f.get_node(f'/{speed_var}').col(f'{speed_val}')
            n_gen = cfg_file['juliet'][s][l]['n_files'] * evts_per_file

            weights = weighting.calc_juliet_weight_from_weight_dict_and_prop_matrix(
                weight_dict=weight_dict, prop_matrix=prop_matrix, 
                flux=gzk_partial, n_gen=n_gen, livetime=livetime
            )

            L2_q_mask = charge > q_cut
            L2_ndoms_mask = ndoms > ndom_cut
            L2_mask = L2_q_mask & L2_ndoms_mask
            
            ehe_weights = np.concatenate((ehe_weights, weights[L2_mask]))
            ehe_charge = np.concatenate((ehe_charge, charge[L2_mask]))
            ehe_speed = np.concatenate((ehe_speed, speed[L2_mask]))
            
            # compute the product of the L2 and L3 mask
            if version=="old":
                L3_mask = analysis_9yr.track_quality_cut_pass_9yr(speed, charge)
            elif version=="new":
                L3_mask = cuts.track_quality_cut(speed, charge)
        
            if do_efficiency:
                # if we want to tabulate efficiencies vs energy
                # then we need to project into enu space

                # no cuts
                the_enu_weights = weighting.calc_juliet_flux_weight(
                    weight_dict = weight_dict, prop_matrix = prop_matrix, 
                    flux= gzk_partial, n_gen = n_gen, livetime=livetime,
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

                # after L3 cuts
                
                L3_mask = np.logical_and(L2_mask, L3_mask)

                the_enu_weights_L3 = weighting.calc_juliet_flux_weight(
                    weight_dict = weight_dict, prop_matrix = prop_matrix, 
                    flux= gzk_partial, n_gen = n_gen, livetime=livetime,
                    selection_mask = L3_mask
                )
                if ehe_energy_weights_after_L3 is None:
                    ehe_energy_weights_after_L3 = copy.deepcopy(the_enu_weights_L3)
                else:
                    ehe_energy_weights_after_L3 += copy.deepcopy(the_enu_weights_L3)

        the_f.close()

if do_efficiency:

    # eff vs energy    
    energy_bins, energy_bin_centers = weighting.get_juliet_enu_binning()

    total_ehe_weight_nocuts = ehe_energy_weights_nocuts.sum()
    total_ehe_weight_L2 = ehe_energy_weights_after_L2.sum()
    total_ehe_weight_L3 = ehe_energy_weights_after_L3.sum()
    total_eff_L2 = total_ehe_weight_L2/total_ehe_weight_nocuts
    total_eff_L3 = total_ehe_weight_L3/total_ehe_weight_nocuts

    fig, ax = plt.subplots(1, 1, figsize=(5,5))

    # eff vs energy
    ax.plot(energy_bin_centers, ehe_energy_weights_after_L2/ehe_energy_weights_nocuts,
        linewidth=3, label="After L2, {:.2f}% Eff".format(total_eff_L2*100.))
    ax.plot(energy_bin_centers, ehe_energy_weights_after_L3/ehe_energy_weights_nocuts,
        linewidth=3, ls="--", label="After L3, {:.2f}% Eff".format(total_eff_L3*100.))
    
    print("Total L2 Eff {:.2f}".format(total_eff_L2*100.))
    print("Total L3 Eff {:.2f}".format(total_eff_L3*100.))
    ax.set_xscale('log')
    ax.set_xlabel(r"E$_{\nu}$ / GeV")
    ax.set_ylabel("Efficiency")
    ax.set_ylim([0, 1.1])
    ax.legend()

    fig.tight_layout()
    fig.savefig("plots/eff_L2_L3_{}.png".format(version), dpi=300)
    del fig, ax


#############################
# 2D histogram of charge and speed, 
# to see how our pivot point works out
#############################

if do_plots:
        
    def get_track_quality_cut_9yr(fit_quality):
        if fit_quality <= 80:
            return 10**4.6
        elif fit_quality > 80 and fit_quality < 120:
            return np.power(10, 4.6 + (0.6/40) * (fit_quality - 80))
        elif fit_quality >= 120:
            return 10**5.2
    
    def get_track_quality_cut_NextGen(speed):
        if speed < 0.26:
            return 10**5.25
        elif speed >= 0.26 and speed < 0.28:
            # return np.power(10, 4.65 + (0.6/0.02) * (speed - 0.26)) 
            return np.power(10, 5.2 - (0.6/0.02) * (speed - 0.26))
        elif speed >= 0.28:
            return 10**4.65
    
    charge_cut_vals = []    
    for s in speed_bins:
        if version == 'old':
            charge_cut_vals.append(get_track_quality_cut_9yr(s))
        elif version == 'new':
            charge_cut_vals.append(get_track_quality_cut_NextGen(s))

    if version == "old":
        pass_mask_cor = analysis_9yr.track_quality_cut_pass_9yr(cor_speed, cor_charge)
        pass_mask_nugen = analysis_9yr.track_quality_cut_pass_9yr(nugen_speed, nugen_charge)
        pass_mask_ehe = analysis_9yr.track_quality_cut_pass_9yr(ehe_speed, ehe_charge)
        print("L3 Pass rate corsika {}".format(np.sum(cor_weights[pass_mask_cor])))
        print("L3 Pass rate nugen {}".format(np.sum(nugen_atmo_weights[pass_mask_nugen])))
        print("L3 Pass rate ehe {}".format(np.sum(ehe_weights[pass_mask_ehe])))
    if version == "new":
        pass_mask_cor = cuts.track_quality_cut(cor_speed, cor_charge)
        pass_mask_nugen = cuts.track_quality_cut(nugen_speed, nugen_charge)
        pass_mask_ehe = cuts.track_quality_cut(ehe_speed, ehe_charge)
        print("L3 Pass rate corsika {}".format(np.sum(cor_weights[pass_mask_cor])))
        print("L3 Pass rate nugen {}".format(np.sum(nugen_atmo_weights[pass_mask_nugen])))
        print("L3 Pass rate ehe {}".format(np.sum(ehe_weights[pass_mask_ehe])))

    # energy distrbution of nugen events
    fig, ax = plt.subplots(1, 1, figsize=(7,5))
    bins = np.logspace(4, 9, 21)
    ax.hist(nugen_energy, bins=bins, weights=nugen_atmo_weights,
        linewidth=3, histtype='step', 
        label="After L2, {:.1e} evts".format(np.sum(nugen_atmo_weights))
    )
    ax.hist(nugen_energy[pass_mask_nugen], bins=bins, weights=nugen_atmo_weights[pass_mask_nugen],
        linewidth=3, histtype='step', linestyle='--',
        label="After L3, {:.1e} evts".format(np.sum(nugen_atmo_weights[pass_mask_nugen]))
    )
    ax.set_ylabel("Weighted Atm Nu Events, {:.2f} days".format(livetime/60/60/24))
    ax.set_xlabel("Neutrino Energy")
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_ylim([1E-6, 1E-1])
    ax.legend()
    fig.savefig("num_events_atmo_nu_L3_{}.png".format(version), dpi=300)
    del fig, ax
        
    plotting_options = {
        'xlabel': r'LineFit Speed',
        'ylabel': 'Q / PE',
        'zlabel': 'Number of Events',
        'cmap': plt.cm.plasma ,
        'norm': colors.LogNorm(),
        'zlims':  (1E-5, 1E2)
    }

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18,5))
    
    # corsika
    cor_sum, cor_xedges, cor_yedges, cor_im = ax1.hist2d(
        x = cor_speed, y = cor_charge, 
        weights = cor_weights, bins = [speed_bins, charge_bins],
        cmap = plotting_options['cmap'], norm=colors.LogNorm()
    )
    ax1.set_title("Atm Muon (Corsika)")
    cor_cbar = plt.colorbar(cor_im, ax=ax1)
    
    # nugen
    nugen_sum, nugen_xedges, nugen_yedges, nugen_im = ax2.hist2d(
        x = nugen_speed, y = nugen_charge, 
        weights = nugen_atmo_weights, bins = [speed_bins, charge_bins],
        cmap = plotting_options['cmap'], norm=colors.LogNorm()
    )
    ax2.set_title("Atm Nu (Nugen)")
    nugen_cbar = plt.colorbar(nugen_im, ax=ax2)
    
    # juliet
    ehe_sum, ehe_xedges, ehe_yedges, ehe_im = ax3.hist2d(
        x = ehe_speed, y = ehe_charge, 
        weights = ehe_weights, bins = [speed_bins, charge_bins],
        cmap = plotting_options['cmap'], norm=colors.LogNorm()
    )
    ax3.set_title("Cosmogenic (Juliet)")
    ehe_cbar = plt.colorbar(ehe_im, ax=ax3)
    
    # decorate all the stuff with titles, blah blah
    for a in [ax1, ax2, ax3]: # axes
        a.set_yscale('log')
        a.set_xlabel(plotting_options['xlabel'])
        a.set_ylabel(plotting_options['ylabel'])
        a.set_xlim(speed_lims)
        a.plot(speed_bins, charge_cut_vals, color='red')
    for i in [cor_im, nugen_im, ehe_im]: # im objects
        i.set_clim(*plotting_options['zlims'])
    for c in [cor_cbar, nugen_cbar, ehe_cbar]: # cbars
        c.set_label(plotting_options['zlabel'])
        
    fig.tight_layout()
    fig.savefig('hist2d_q_speed_{}.png'.format(version), dpi=300)

    del fig, ax1, ax2, ax3
