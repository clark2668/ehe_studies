import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import style
import tables
import pandas as pd
import copy
import yaml
import numpy as np

from functools import partial
from ehefluxes import fluxes
from eheanalysis import weighting, plotting
gzk_flux_18 = fluxes.EHEFlux("cosmogenic_ahlers2010_1E18")
gzk_flux_19 = fluxes.EHEFlux("cosmogenic_ahlers2010_1E19")
gzk_partial_18 = partial(gzk_flux_18, which_species="nue_sum") 
gzk_partial_19 = partial(gzk_flux_19, which_species="nue_sum") 
cr_flux = weighting.get_flux_model('GaisserH4a', 'corsika')
atmo_flux = weighting.get_flux_model('H3a_SIBYLL23C', 'nugen')

def astro_flux(energy):
    # flux of mu @ 100 TeV (per-flavor, per-particle flux)
    # then multiply by two to get the nu+nubar sum
    return 2 * 1.44e-18 / 2 * (energy/1e5)**-2.37


style.use('/home/brian/IceCube/ehe/ehe_software/ehe_code/EHE_analysis/eheanalysis/ehe.mplstyle')

livetime = 365 * 24 * 60 * 60
print(livetime)

cfg_file = 'config.yaml'
cfg_file = yaml.safe_load(open(cfg_file))

which_method='old'
if which_method == 'new':
    charge_var = cfg_file['variables']['hqtot']['variable']
    charge_val = cfg_file['variables']['hqtot']['value']
    ndoms_var =  cfg_file['variables']['ndoms']['variable']
    ndoms_val =  cfg_file['variables']['ndoms']['value']
    classifier_var = cfg_file['variables']['speed']['variable']
    classifier_val = cfg_file['variables']['speed']['value']
    classifier_bins = np.linspace(0, 2, 601)
    classifier_label = "LineFit Speed"
    classifier_cut = 0.27
    classifier_lims = [0, 0.5]
    qcut = 27500
    log10_q_cut = np.log10(qcut)
    nue_cumulative_sign = -1
elif which_method == 'old':
    charge_var = "EHEPortiaEventSummarySRT"
    charge_val = "bestNPEbtw"
    ndoms_var =  "EHEPortiaEventSummarySRT"
    ndoms_val =  "NCHbtw"
    classifier_var = "EHEOpheliaSRT_ImpLF"
    classifier_val = "fitQuality"
    classifier_label = "Ophelia FitQual"
    classifier_bins = np.linspace(0, 1000, 601)
    classifier_cut = 100
    classifier_lims = [0,500]
    qcut = 25000
    log10_q_cut = np.log10(qcut)
    nue_cumulative_sign = 1
ndom_cut = 100


ehe_weights_gzk18 ={
    "nue": np.asarray([]),
    "mu": np.asarray([])
}
ehe_weights_gzk19 ={
    "nue": np.asarray([]),
    "mu": np.asarray([])
}
ehe_weights_astro ={
    "nue": np.asarray([]),
    "mu": np.asarray([])
}

ehe_charge = {
    "nue": np.asarray([]),
    "mu": np.asarray([])
}
ehe_ndoms = {
    "nue": np.asarray([]),
    "mu": np.asarray([])
}
ehe_classifier = {
    "nue": np.asarray([]),
    "mu": np.asarray([])
}


juliet_species = ["nue", "mu"]
juliet_energy_levels = ["high_energy", "very_high_energy"]
# juliet_energy_levels = ["high_energy"]
# juliet_energy_levels = []
events_per_file = {
    "nue_high_energy": 600,
    "nue_very_high_energy": 80,
    "mu_high_energy": 150,
    "mu_very_high_energy": 20
}

#############################
# ehe/cosmogenic flux (juliet)
#############################


for s in juliet_species:
    for l in juliet_energy_levels:

        print(f"Working on juliet {s} {l}")

        the_f = tables.open_file(cfg_file['juliet'][s][l]['file'])
        weight_dict, prop_matrix, evts_per_file = weighting.get_juliet_weightdict_and_propmatrix(the_f)

        evts_per_file = events_per_file[f"{s}_{l}"] # override to fix L2 issue

        charge = the_f.get_node(f"/{charge_var}").col(f"{charge_val}")
        classifier = the_f.get_node(f"/{classifier_var}").col(f"{classifier_val}")
        n_gen = cfg_file['juliet'][s][l]['n_files'] * evts_per_file
        ndoms = the_f.get_node(f'/{ndoms_var}').col(f'{ndoms_val}')

        weight_gzk18 = weighting.calc_juliet_weight_from_weight_dict_and_prop_matrix(
            weight_dict=weight_dict, prop_matrix=prop_matrix, 
            flux=gzk_partial_18, 
            n_gen=n_gen, livetime=livetime,
        )
        # weight_gzk19 = weighting.calc_juliet_weight_from_weight_dict_and_prop_matrix(
        #     weight_dict=weight_dict, prop_matrix=prop_matrix, 
        #     flux=gzk_partial_19, 
        #     n_gen=n_gen, livetime=livetime,
        # )
        # weight_astro = weighting.calc_juliet_weight_from_weight_dict_and_prop_matrix(
        #     weight_dict=weight_dict, prop_matrix=prop_matrix, 
        #     flux=astro_flux, 
        #     n_gen=n_gen, livetime=livetime,
        # )        

        ehe_weights_gzk18[s] = np.concatenate((ehe_weights_gzk18[s], copy.deepcopy(abs(weight_gzk18))))
        # ehe_weights_gzk19[s] = np.concatenate((ehe_weights_gzk19[s], copy.deepcopy(abs(weight_gzk19))))
        # ehe_weights_astro[s] = np.concatenate((ehe_weights_astro[s], copy.deepcopy(abs(weight_astro))))
        ehe_charge[s] = np.concatenate((ehe_charge[s], copy.deepcopy(charge)))
        ehe_classifier[s] = np.concatenate((ehe_classifier[s], copy.deepcopy(classifier)))
        ehe_ndoms[s] = np.concatenate((ehe_ndoms[s], copy.deepcopy(ndoms)))
        the_f.close()

ehe_mask_L1 = {
    "nue": ehe_charge['nue']>0,
    "mu": ehe_charge['mu']>0
}

ehe_mask_L2 = {
    "nue": np.logical_and(ehe_charge['nue']>qcut, ehe_ndoms['nue']>ndom_cut),
    "mu":  np.logical_and(ehe_charge['mu']>qcut,  ehe_ndoms['mu']>ndom_cut)
}


#############################
# muon bundles (corsika)
#############################

cor_classifier = np.asarray([])
cor_charge = np.asarray([])
cor_weights = np.asarray([])
cor_ndoms = np.asarray([])

corsika_sets = ["20787"]
for c in corsika_sets:
    print("Working on corsika {}".format(c))
    cor_file = pd.HDFStore(cfg_file['corsika'][c]['file'])
    cor_weighter = weighting.get_weighter( cor_file, 'corsika',
        cfg_file['corsika'][c]['n_files']
        )
    this_cor_charge = cor_weighter.get_column(charge_var, charge_val)
    this_cor_classifier = cor_weighter.get_column(classifier_var, classifier_val)
    this_cor_ndoms = cor_weighter.get_column(ndoms_var, ndoms_val)
    
    cor_charge = np.concatenate((cor_charge, this_cor_charge))
    cor_classifier = np.concatenate((cor_classifier, this_cor_classifier))
    cor_ndoms = np.concatenate((cor_ndoms, this_cor_ndoms))
    cor_weights = np.concatenate((cor_weights, cor_weighter.get_weights(cr_flux) * livetime))
    cor_file.close()


cor_mask_L1 = cor_charge > 0
cor_mask_L2 = np.logical_and( cor_charge > qcut, cor_ndoms > ndom_cut)

if which_method is 'new':
    cor_track_mask = cor_classifier >= classifier_cut
    # need to build the cascade mask first, since we overwrite the track def
    cor_cascade_mask = np.logical_and(~cor_track_mask, cor_mask_L2)
    cor_track_mask = np.logical_and(cor_track_mask, cor_mask_L2)
elif which_method is 'old':
    cor_track_mask = cor_classifier <= classifier_cut
    cor_cascade_mask = np.logical_and(~cor_track_mask, cor_mask_L2)
    cor_track_mask = np.logical_and(cor_track_mask, cor_mask_L2)

cor_weight_total = np.sum(cor_weights[cor_mask_L2])
cor_weight_tracks = np.sum(cor_weights[cor_track_mask])
cor_weight_cascades = np.sum(cor_weights[cor_cascade_mask])
cor_misclassifier = cor_weight_cascades/cor_weight_total
print(f"Cor Correct: {cor_weight_tracks/cor_weight_total:.3f}")
print(f"Cor False: {cor_misclassifier:.3f}")

def scaling(ax):
    ax.set_xlim(classifier_lims)
    ax.axvline(classifier_cut, linestyle='--')
lw = 2


make_1d_plots_diffcutlevels = False
if make_1d_plots_diffcutlevels:
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10,5))
    
    n, b, p = ax1.hist( 
        ehe_classifier['mu'][ehe_mask_L1['mu']], bins=classifier_bins, 
        weights=ehe_weights_gzk18["mu"][ehe_mask_L1['mu']],
        histtype='step', label=r'GZK $\mu$, L1', linewidth=lw, density=True)
    n, b, p = ax1.hist( 
        ehe_classifier['mu'][ehe_mask_L2['mu']], bins=classifier_bins, 
        weights=ehe_weights_gzk18["mu"][ehe_mask_L2['mu']],
        histtype='step', label=r'GZK $\mu$, L2', linewidth=lw, density=True)
    ax1.set_yscale('log')
    ax1.set_xlabel(classifier_label)
    ax1.set_ylabel('Events / {:.2f} days'.format(livetime/(60*60*24)))
    ax1.set_ylim([1E-5, 1E2])
    ax1.legend(loc='upper left')
    
    n, b, p = ax2.hist( 
        cor_classifier[cor_mask_L1], bins=classifier_bins, 
        weights=cor_weights[cor_mask_L1],
        histtype='step', label=r'Corsika, L1', linewidth=lw, density=True)
    n, b, p = ax2.hist( 
        cor_classifier[cor_mask_L2], bins=classifier_bins, 
        weights=cor_weights[cor_mask_L2],
        histtype='step', label=r'Corsika, L2', linewidth=lw, density=True)
    ax2.set_yscale('log')
    ax2.set_xlabel(classifier_label)
    ax2.set_ylabel('Events / {:.2f} days'.format(livetime/(60*60*24)))
    ax2.set_ylim([1E-5, 1E2])
    ax2.legend(loc='upper left')    
    
    scaling(ax1)
    scaling(ax2)

    fig.tight_layout()
    fig.savefig(f'./figs/classifier_{which_method}_L1vsL2.png', dpi=300)

# need to also plot speed vs E_primary (from Polyplopia Primary)
# need to provde that it's these high energy events that are cutting in
# and creating the bump, and they're simply not
# in the corsika


make_1d_plots_difffluxes = False
if make_1d_plots_difffluxes:
    
    # plot all L1 together, plot all L2 together
    fig, ((ax3, ax4), (ax1, ax2)) = plt.subplots(2, 2, figsize=(10,10))
    
    # L1 (two fluxes together)
    n, b, p = ax1.hist( 
        ehe_classifier['mu'][ehe_mask_L1['mu']], bins=classifier_bins, 
        weights=ehe_weights_gzk18["mu"][ehe_mask_L1['mu']],
        histtype='step', label=r'GZK E18 $\mu$, L1', linewidth=lw, density=True)
    n, b, p = ax1.hist( 
        ehe_classifier['mu'][ehe_mask_L1['mu']], bins=classifier_bins, 
        weights=ehe_weights_astro["mu"][ehe_mask_L1['mu']],
        histtype='step', label=r'Astro $\mu$, L1', linewidth=lw, density=True)
    
    # L2 (two fluxes together)
    n, b, p = ax2.hist( 
        ehe_classifier['mu'][ehe_mask_L2['mu']], bins=classifier_bins, 
        weights=ehe_weights_gzk18["mu"][ehe_mask_L2['mu']],
        histtype='step', label=r'GZK E18 $\mu$, L2', linewidth=lw, density=True)
    n, b, p = ax2.hist( 
        ehe_classifier['mu'][ehe_mask_L2['mu']], bins=classifier_bins, 
        weights=ehe_weights_astro["mu"][ehe_mask_L2['mu']],
        histtype='step', label=r'Astro $\mu$, L2', linewidth=lw, density=True)

    # L1 & L2 for GZK together
    n, b, p = ax3.hist( 
        ehe_classifier['mu'][ehe_mask_L1['mu']], bins=classifier_bins, 
        weights=ehe_weights_gzk18["mu"][ehe_mask_L1['mu']],
        histtype='step', label=r'GZK E18 $\mu$, L1', linewidth=lw, density=True)
    n, b, p = ax3.hist( 
        ehe_classifier['mu'][ehe_mask_L2['mu']], bins=classifier_bins, 
        weights=ehe_weights_gzk18["mu"][ehe_mask_L2['mu']],
        histtype='step', label=r'GZK E18 $\mu$, L2', linewidth=lw, density=True)
    
    # L1 & L2 for astro together
    n, b, p = ax4.hist( 
        ehe_classifier['mu'][ehe_mask_L1['mu']], bins=classifier_bins, 
        weights=ehe_weights_astro["mu"][ehe_mask_L1['mu']],
        histtype='step', label=r'Astro $\mu$, L1', linewidth=lw, density=True)
    n, b, p = ax4.hist( 
        ehe_classifier['mu'][ehe_mask_L2['mu']], bins=classifier_bins, 
        weights=ehe_weights_astro["mu"][ehe_mask_L2['mu']],
        histtype='step', label=r'Astro $\mu$, L2', linewidth=lw, density=True)

    for ax in (ax1, ax2, ax3, ax4):
        ax.set_yscale('log')
        ax.set_xlabel(classifier_label)
        ax.set_ylabel('PDF')
        ax.set_ylim([1E-2, 1E2])
        ax.legend(loc='upper left')
        scaling(ax)
    
    fig.tight_layout()
    fig.savefig(f'./figs/classifier_{which_method}_L1vsL2_difffluxes.png', dpi=300)


make_1d_plots_3panel = False
if make_1d_plots_3panel:
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15,5))

    ax1.hist( 
        ehe_classifier['nue'][ehe_mask_L2['nue']], bins=classifier_bins, 
        weights=ehe_weights_gzk18["nue"][ehe_mask_L2['nue']],
        histtype='step', label=r'GZK $\nu_{e}$', linewidth=lw)
    n, b, p = ax1.hist( 
        ehe_classifier['mu'][ehe_mask_L2['mu']], bins=classifier_bins, 
        weights=ehe_weights_gzk18["mu"][ehe_mask_L2['mu']],
        histtype='step', label=r'GZK $\mu$', linewidth=lw)
    ax1.hist( cor_classifier[cor_mask_L2], bins=classifier_bins, weights=cor_weights[cor_mask_L2],
            histtype='step', label='Corsika (H4a)', linewidth=lw)
    ax1.set_yscale('log')
    ax1.set_xlabel(classifier_label)
    ax1.set_ylabel('Events / {:.2f} days'.format(livetime/(60*60*24)))
    ax1.set_ylim([1E-6, 1E4])
    ax1.legend()

    ax2.hist( 
        ehe_classifier['nue'][ehe_mask_L2['nue']], bins=classifier_bins, 
        weights=ehe_weights_gzk18["nue"][ehe_mask_L2['nue']],
        histtype='step', linewidth=lw, density=True)
    ax2.hist( 
        ehe_classifier['mu'][ehe_mask_L2['mu']], bins=classifier_bins, 
        weights=ehe_weights_gzk18["mu"][ehe_mask_L2['mu']],
        histtype='step', linewidth=lw, density=True)
    ax2.hist( 
        cor_classifier[cor_mask_L2], bins=classifier_bins, 
        weights=cor_weights[cor_mask_L2],
        histtype='step', linewidth=lw, density=True)

    ax2.set_xlabel(classifier_label)
    ax2.set_ylabel('Density')

    ax3.hist( 
        ehe_classifier['nue'][ehe_mask_L2['nue']], bins=classifier_bins, 
        weights=ehe_weights_gzk18["nue"][ehe_mask_L2['nue']],
        histtype='step', linewidth=lw, density=True, cumulative=nue_cumulative_sign)
    ax3.hist( 
        ehe_classifier['mu'][ehe_mask_L2['mu']], bins=classifier_bins, 
        weights=ehe_weights_gzk18["mu"][ehe_mask_L2['mu']],
        histtype='step', linewidth=lw, density=True, cumulative=-nue_cumulative_sign)
    ax3.hist( 
        cor_classifier[cor_mask_L2], bins=classifier_bins, 
        weights=cor_weights[cor_mask_L2],
        histtype='step', linewidth=lw, density=True, cumulative=-nue_cumulative_sign)
    ax3.set_xlabel(classifier_label)
    ax3.set_ylabel('CDF')

    expected_rate = np.sum(ehe_weights_gzk18['mu'][ehe_mask_L2['mu']])
    print(f"Expected mu rate: {expected_rate}")

    for ax in [ax1, ax2, ax3]:
        scaling(ax)

    fig.tight_layout()
    fig.savefig(f'./figs/classifier_{which_method}.png', dpi=300)


make_confusion = True
if make_confusion:
    def score_classifier(classifier_results, version='new'):
        # 0 means cascade, 1 means track
        copy = classifier_results.copy()
        if version is 'new':
            copy[copy<0.27] =  0 # cascade
            copy[copy>=0.27] = 1 # track
        elif version is 'old':
            copy[copy<=100] =  1 # track
            copy[copy>100] =   0 # cascade
        return copy

    scores = {
        "nue": score_classifier(ehe_classifier['nue'], which_method),
        "mu": score_classifier(ehe_classifier['mu'], which_method)
    }

    y_pred = np.concatenate((
        scores['nue'][ehe_mask_L2['nue']], 
        scores['mu'][ehe_mask_L2['mu']]
        ))
    y_true = np.concatenate((
        np.full_like(ehe_classifier['nue'], 0)[ehe_mask_L2['nue']],
        np.full_like(ehe_classifier['mu'],  1)[ehe_mask_L2['mu']]
    ))
    y_weights = np.concatenate((
        ehe_weights_gzk18['nue'][ehe_mask_L2['nue']],
        ehe_weights_gzk18['mu'][ehe_mask_L2['mu']]
    ))

    from sklearn import metrics
    confusion_matrix = metrics.confusion_matrix(
        y_true, y_pred, 
        sample_weight=y_weights, 
        normalize='true'
        )
    cm_display = metrics.ConfusionMatrixDisplay(
        confusion_matrix=confusion_matrix,
        display_labels=['Cascade', 'Track']
        )
    cm_display.plot()
    plt.tight_layout()
    plt.savefig(f'./figs/confusion_{which_method}.png', dpi=300)
