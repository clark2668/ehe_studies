'''
Calculate the histogram of passing neutrinos vs energy
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
version = 'new'
version = 'old'
if version == 'new':
    charge_var = cfg_file['variables']['charge']['variable']
    charge_val = cfg_file['variables']['charge']['value']
    ndoms_var =  cfg_file['variables']['ndoms']['variable']
    ndoms_val =  cfg_file['variables']['ndoms']['value']
    speed_var = cfg_file['variables']['speed']['variable']
    speed_val = cfg_file['variables']['speed']['value']
    zen_var = cfg_file['variables']['zenith']['variable']
    zen_val = cfg_file['variables']['zenith']['value']
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
    zen_var = "EHEOpheliaParticleSRT_ImpLF"
    zen_val = "zenith"
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

gzk_flux_1EeV = fluxes.EHEFlux("cosmogenic_ahlers2010_1E18")
gzk_partial_1EeV = partial(gzk_flux_1EeV, which_species="nue_sum") 
gzk_flux_3EeV = fluxes.EHEFlux("cosmogenic_ahlers2010_1E185")
gzk_partial_3EeV = partial(gzk_flux_3EeV, which_species="nue_sum") 
gzk_flux_10EeV = fluxes.EHEFlux("cosmogenic_ahlers2010_1E19")
gzk_partial_10EeV = partial(gzk_flux_10EeV, which_species="nue_sum") 

style.use('/home/brian/IceCube/ehe/max_tools/EHE_analysis/eheanalysis/ehe.mplstyle')

livetime = 3142.5 * 24 * 60 * 60 # number of seconds in the 9yr paper livetime

juliet_species = ["nue", "numu", "nutau", "mu", "tau"]
juliet_energy_levels = ["high_energy", "very_high_energy"]

ehe_energy_weights_1EeV = None
ehe_energy_weights_3EeV = None
ehe_energy_weights_10EeV = None

for s in juliet_species:
    for l in juliet_energy_levels:

        if s in juliet_species:

            print(f"Working on juliet {s}, {l}")
            the_f = tables.open_file(cfg_file['juliet'][s][l]['file'])
            weight_dict, prop_matrix, evts_per_file = weighting.get_juliet_weightdict_and_propmatrix(the_f)
            charge = the_f.get_node(f'/{charge_var}').col(f'{charge_val}')
            ndoms = the_f.get_node(f'/{ndoms_var}').col(f'{ndoms_val}')
            speed = the_f.get_node(f'/{speed_var}').col(f'{speed_val}')
            zen = the_f.get_node(f'/{zen_var}').col(f'{zen_val}')
            n_gen = cfg_file['juliet'][s][l]['n_files'] * evts_per_file

            L2_q_mask = charge > q_cut
            L2_ndoms_mask = ndoms > ndom_cut
            L2_mask = L2_q_mask & L2_ndoms_mask
                       
            # compute the product of the L2 and L3 mask
            if version=="old":
                L3_mask = analysis_9yr.track_quality_cut_pass_9yr(speed, charge)
                L4_mask = analysis_9yr.muon_bundle_cut_pass_9yr(zen, charge)
            elif version=="new":
                L3_mask = cuts.track_quality_cut(speed, charge)
                L4_mask = cuts.track_quality_cut(speed, charge)
            
            L3_mask = np.logical_and(L2_mask, L3_mask)
            L4_mask = np.logical_and(L3_mask, L4_mask)

            the_enu_weights_1EeV = weighting.calc_juliet_flux_weight(
                weight_dict = weight_dict, prop_matrix = prop_matrix, 
                flux= gzk_partial_1EeV, n_gen = n_gen, livetime=livetime,
                selection_mask = L4_mask
            )
            the_enu_weights_3EeV = weighting.calc_juliet_flux_weight(
                weight_dict = weight_dict, prop_matrix = prop_matrix, 
                flux= gzk_partial_3EeV, n_gen = n_gen, livetime=livetime,
                selection_mask = L4_mask
            )
            the_enu_weights_10EeV = weighting.calc_juliet_flux_weight(
                weight_dict = weight_dict, prop_matrix = prop_matrix, 
                flux= gzk_partial_10EeV, n_gen = n_gen, livetime=livetime,
                selection_mask = L4_mask
            )

            if ehe_energy_weights_1EeV is None:
                ehe_energy_weights_1EeV = copy.deepcopy(the_enu_weights_1EeV)
                ehe_energy_weights_3EeV = copy.deepcopy(the_enu_weights_3EeV)
                ehe_energy_weights_10EeV = copy.deepcopy(the_enu_weights_10EeV)

            else:
                ehe_energy_weights_1EeV += copy.deepcopy(the_enu_weights_1EeV)
                ehe_energy_weights_3EeV += copy.deepcopy(the_enu_weights_3EeV)
                ehe_energy_weights_10EeV += copy.deepcopy(the_enu_weights_10EeV)
        
        the_f.close()

energy_bins, energy_bin_centers = weighting.get_juliet_enu_binning()

fig, ax =  plt.subplots(1,1, figsize=(8,5))
n, b, p = ax.hist(energy_bin_centers, bins=energy_bins, weights=ehe_energy_weights_1EeV,
                  histtype='step', label=f'Total = {ehe_energy_weights_1EeV.sum():.2f}, 1 EeV',
                  linewidth=3)
n, b, p = ax.hist(energy_bin_centers, bins=energy_bins, weights=ehe_energy_weights_3EeV,
                  histtype='step', label=f'Total = {ehe_energy_weights_3EeV.sum():.2f}, 3 EeV',
                  linewidth=3)
n, b, p = ax.hist(energy_bin_centers, bins=energy_bins, weights=ehe_energy_weights_10EeV,
                  histtype='step', label=f'Total = {ehe_energy_weights_10EeV.sum():.2f}, 10 EeV',
                  linewidth=3)
ax.legend()
ax.set_ylabel("Events / {:.2f} days".format(livetime/60/60/24))
ax.set_xlabel(r'Energy [GeV]')
ax.set_ylabel('Counts')
ax.set_xscale('log')
fig.tight_layout()
fig.savefig(f'hist_passing.png', dpi=300)
del fig, ax


