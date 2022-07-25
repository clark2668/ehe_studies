'''
Point of this script is to understand how well our proposed
LineFit variable works in Juliet
'''

import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import style
import tables
import pandas as pd
import copy
import yaml

from eheanalysis import plotting

cfg_file = 'config.yaml'
cfg_file = yaml.safe_load(open(cfg_file))

alt_ndoms_var = 'HitMultiplicityValues'


version = 'new'
if version == 'new':
    charge_var = cfg_file['variables']['charge']['variable']
    charge_val = cfg_file['variables']['charge']['value']
    ndoms_var =  cfg_file['variables']['ndoms']['variable']
    ndoms_val =  cfg_file['variables']['ndoms']['value']
    speed_var = cfg_file['variables']['speed']['variable']
    speed_val = cfg_file['variables']['speed']['value']
    speed_bins = np.linspace(0, 2, 201)
    speed_label = "LineFit Speed"
    speed_cut = 0.26
    speed_lims = [0, 0.5]
    log10_q_cut = np.log10(27500)
elif version == 'old':
    charge_var = "EHEPortiaEventSummarySRT"
    charge_val = "bestNPEbtw"
    ndoms_var =  "EHEPortiaEventSummarySRT"
    ndoms_val =  "NCHbtw"
    speed_var = "EHEOpheliaSRT_ImpLF"
    speed_val = "fitQuality"
    speed_label = "Ophelia FitQual"
    speed_bins = np.linspace(0, 600, 601)
    speed_cut = 100
    speed_lims = [0,500]
    log10_q_cut = np.log10(25000)

speed_bin_centers = plotting.get_bin_centers(speed_bins)
q_cut = np.power(10., log10_q_cut)
ndom_cut = 100   


from functools import partial
from ehefluxes import fluxes
from eheanalysis import weighting
gzk_flux = fluxes.EHEFlux("ahlers_gzk")
gzk_partial = partial(gzk_flux, which_species="nue_sum") 
cr_flux = weighting.get_flux_model('GaisserH4a', 'corsika')

# set up datasets
# juliet (EHE neutrinos)
juliet_species = ["nue", "numu", "nutau", "mu", "tau"]
juliet_energy_levels = ["high_energy"]
juliet_species = ["nue"]

corsika_sets = ["20787"]

burn_samples = ["IC86-I-pass2", "IC86-II-pass2", "IC86-III-pass2"]

style.use('/home/brian/IceCube/ehe/max_tools/EHE_analysis/eheanalysis/ehe.mplstyle')

livetime = 0
for b in burn_samples:
    livetime += cfg_file['burn_sample'][b]['livetime']
print("Total livetime {}".format(livetime))



#############################
# ehe/cosmogenic flux (juliet)
#############################

# we want EHE speed broken out a little more carefully...
# maybe nue, numu, corsika?

summed_ehe_nue = None
summed_ehe_mu = None

juliet_species = ["nue", "mu"]
for s in juliet_species:
    for l in juliet_energy_levels:

        if s in juliet_species:

            print(f"Working on juliet {s}")
            the_f = tables.open_file(cfg_file['juliet'][s][l]['file'])
            weight_dict, prop_matrix, evts_per_file = weighting.get_juliet_weightdict_and_propmatrix(the_f)
            charge = the_f.get_node(f'/{charge_var}').col(f'{charge_val}')
            ndoms = the_f.get_node(f'/{ndoms_var}').col(f'{ndoms_val}')
            speed = the_f.get_node(f'/{speed_var}').col(f'{speed_val}')
            n_gen = cfg_file['juliet'][s][l]['n_files'] * evts_per_file

            L2_q_mask = charge > q_cut
            L2_ndom_mask = ndoms > ndom_cut
            L2_mask = L2_q_mask & L2_ndom_mask # build the L2 cut

            h_speed = weighting.make_enu_2d_hist( weight_dict, n_gen, gzk_partial,
                var_values = speed, var_bins = speed_bins,
                prop_matrix = prop_matrix, livetime=livetime,
                selection_mask = L2_mask
            )
            
            if s in ["nue"]:
                if summed_ehe_nue is None:
                    summed_ehe_nue = copy.deepcopy(h_speed)
                else:
                    summed_ehe_nue += copy.deepcopy(h_speed)
            else:
                if summed_ehe_mu  is None:
                    summed_ehe_mu = copy.deepcopy(h_speed)
                else:
                    summed_ehe_mu += copy.deepcopy(h_speed)        
        the_f.close()

ehe_speed_nue = np.asarray([])
ehe_speed_mu = np.asarray([])

if summed_ehe_nue is not None:
    ehe_speed_nue = summed_ehe_nue.sum(axis=1) # project down onto the speed axis
if summed_ehe_mu is not None:
    ehe_speed_mu = summed_ehe_mu.sum(axis=1)

#############################
# muon bundles (corsika)
#############################
cor_speed = np.asarray([])
cor_weights = np.asarray([])
for c in corsika_sets:
    print("Working on corsika {}".format(c))
    cor_file = pd.HDFStore(cfg_file['corsika'][c]['file'])
    cor_weighter = weighting.get_weighter( cor_file, 'corsika',
        cfg_file['corsika'][c]['n_files']
        )
    this_cor_charge = cor_weighter.get_column(charge_var, charge_val)
    # special values for corsika (ugh)
    this_cor_ndoms = cor_weighter.get_column('HitMultiplicityValues', ndoms_val)
    this_cor_speed = cor_weighter.get_column('LineFit_redo', speed_val)
    
    L2_q_mask = this_cor_charge > q_cut
    L2_ndom_mask = this_cor_ndoms > ndom_cut
    L2_mask = L2_q_mask & L2_ndom_mask # build the L2 cut

    cor_speed = np.concatenate((cor_speed, this_cor_speed[L2_mask]))
    cor_weights = np.concatenate((cor_weights, cor_weighter.get_weights(cr_flux)[L2_mask] * livetime))
    cor_file.close()

#############################
# histogram of speed
#############################

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, sharex=True, figsize=(15,5))


ax1.hist( speed_bin_centers, bins=speed_bins, weights=ehe_speed_nue,
        histtype='step', label=r'GZK $\nu_{e}$')
ax1.hist( speed_bin_centers, bins=speed_bins, weights=ehe_speed_mu,
        histtype='step', label=r'GZK $\mu$')
ax1.hist( cor_speed, bins=speed_bins, weights=cor_weights,
        histtype='step', label='Corsika (H4a)')
ax1.set_yscale('log')
ax1.set_xlabel(f"{speed_var}")
ax1.set_ylabel('Events / {:.2f} days'.format(livetime/(60*60*24)))
ax1.set_ylim([1E-7, 1E5])
ax1.legend()


n_nue, bins_nue, patches_nue = ax2.hist( speed_bin_centers, bins=speed_bins, weights=ehe_speed_nue,
        histtype='step',  density=True)
ax2.hist( speed_bin_centers, bins=speed_bins, weights=ehe_speed_mu,
        histtype='step', density=True)
n_cor, bins_cor, patches_cor = ax2.hist( cor_speed, bins=speed_bins, weights=cor_weights,
        histtype='step', density=True)
ax2.set_xlabel(f"{speed_var}")
ax2.set_ylabel('Normalized Counts')
ax2.axvline(speed_cut, linstyle='--')

# figure out efficiency
bins_cor = (bins_cor[1:] + bins_cor[:-1])/2
bins_cor_cut = bins_cor > speed_cut
integral = np.sum(n_cor[bins_cor_cut]) * (bins_cor[1]-bins_cor[0])
print("fitspeed: the cor integral is {}".format(integral))

bins_nue = (bins_nue[1:] + bins_nue[:-1])/2
bins_nue_cut = bins_nue < speed_cut
integral = np.sum(n_nue[bins_nue_cut]) * (bins_nue[1]- bins_nue[0])
print("fitspeed: the nue integral is {}".format(integral))


ax3.hist( speed_bin_centers, bins=speed_bins, weights=ehe_speed_nue,
        histtype='step', density=True, cumulative=-1)
ax3.hist( speed_bin_centers, bins=speed_bins, weights=ehe_speed_mu,
        histtype='step', density=True, cumulative=1)
ax3.hist( cor_speed, bins=speed_bins, weights=cor_weights,
        histtype='step', density=True, cumulative=1)
ax3.set_xlabel("LineFit Speed")
ax3.set_ylabel('CDF')
ax3.axvline(speed_cut, linstyle='--')

def scaling(ax):
    ax.set_xlim(speed_lims)
for ax in [ax1, ax2, ax3]:
    scaling(ax)



fig.tight_layout()
fig.savefig(f'speed_hist_{version}.png')
