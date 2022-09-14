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

from eheanalysis import plotting, analysis_9yr

cfg_file = 'config.yaml'
cfg_file = yaml.safe_load(open(cfg_file))

version = 'old'
# version = "new"
if version == 'new':
    charge_var = cfg_file['variables']['charge']['variable']
    charge_val = cfg_file['variables']['charge']['value']
    ndoms_var =  cfg_file['variables']['ndoms']['variable']
    ndoms_val =  cfg_file['variables']['ndoms']['value']
    speed_var = cfg_file['variables']['speed']['variable']
    speed_val = cfg_file['variables']['speed']['value']
    speed_bins = np.linspace(0, 2, 601)
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
    speed_bins = np.linspace(0, 1000, 601)
    speed_cut = 100
    speed_lims = [0,500]
    log10_q_cut = np.log10(25000)
    nue_cumulative_sign = 1

charge_bins = np.logspace(4, 7, 16)
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
atmo_flux = weighting.get_flux_model('H3a_SIBYLL23C_pr', 'nugen')

burn_samples = ["IC86-I-pass2", "IC86-II-pass2", "IC86-III-pass2"]

style.use('/home/brian/IceCube/ehe/max_tools/EHE_analysis/eheanalysis/ehe.mplstyle')

livetime = 0
for b in burn_samples:
    livetime += cfg_file['burn_sample'][b]['livetime']
livetime = 203.8 * 24 * 60 * 60 # seconds (203.8 days to seconds)
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

nugen_sets = ["nue", "numu", "nutau"]
# nugen_sets = ["nue"]
for n in nugen_sets:
    print("Working on nugen {}".format(n))

    nugen_file = pd.HDFStore(cfg_file['nugen'][n]['file'])
    nugen_weighter = weighting.get_weighter( nugen_file,  'nugen',
        cfg_file['nugen'][n]['n_files']
        )
    this_nugen_charge = np.asarray(nugen_weighter.get_column(charge_var, charge_val))    
    this_nugen_ndoms = np.asarray(nugen_weighter.get_column(ndoms_var, ndoms_val))
    this_nugen_speed = nugen_weighter.get_column(speed_var, speed_val)
    
    L2_q_mask = this_nugen_charge > q_cut
    L2_ndom_mask = this_nugen_ndoms > ndom_cut
    L2_mask = L2_q_mask & L2_ndom_mask # build the L2 cut
    
    nugen_charge = np.concatenate((nugen_charge, this_nugen_charge[L2_mask]))       
    nugen_speed = np.concatenate((nugen_speed, this_nugen_speed[L2_mask]))
    
    this_nugen_atmo_weights = nugen_weighter.get_weights(atmo_flux) * livetime
    nugen_atmo_weights = np.concatenate((nugen_atmo_weights, this_nugen_atmo_weights[L2_mask]))

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
juliet_species = ["numu"]
juliet_energy_levels = ["high_energy"]

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
    
        the_f.close()


#############################
# 2D histogram of charge and speed, 
# to see how our pivot point works out
#############################

do_speed_hist_2d = True
if do_speed_hist_2d:
        
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
    ax1.set_title("Atm Nu (Corsika)")
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