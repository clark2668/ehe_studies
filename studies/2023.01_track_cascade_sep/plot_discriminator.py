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
gzk_flux = fluxes.EHEFlux("cosmogenic_ahlers2010_1E18")
gzk_partial = partial(gzk_flux, which_species="nue_sum") 
cr_flux = weighting.get_flux_model('GaisserH4a', 'corsika')
atmo_flux = weighting.get_flux_model('H3a_SIBYLL23C', 'nugen')

style.use('/home/brian/IceCube/ehe/ehe_software/ehe_code/EHE_analysis/eheanalysis/ehe.mplstyle')

livetime = 365 * 24 * 60 * 60
print(livetime)

cfg_file = 'config.yaml'
cfg_file = yaml.safe_load(open(cfg_file))

which_method='new'
if which_method is 'new':
    qcut = 27500
    charge_var = 'hqtot'
    classifier_var = 'speed'
    classifier_bins = np.linspace(0, 0.5, 150)
    nue_cumulative_sign = -1
    classifier_lims = [0, 0.5]
    classifier_cut = 0.27
if which_method is 'old':
    qcut = 25000
    charge_var = 'portia'
    classifier_var = 'redchisqu'
    classifier_bins = np.linspace(0, 1000, 601)
    nue_cumulative_sign = 1 
    classifier_lims = [0, 500]
    classifier_cut = 100



ehe_weights ={
    "nue": np.asarray([]),
    "mu": np.asarray([])
}
ehe_charge = {
    "nue": np.asarray([]),
    "mu": np.asarray([])
}
ehe_classifier = {
    "nue": np.asarray([]),
    "mu": np.asarray([])
}


juliet_species = ["nue", "mu"]
juliet_energy_levels = ["high_energy", "very_high_energy"]
juliet_energy_levels = ["high_energy"]
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

        charge = the_f.get_node(f"/{cfg_file['variables'][charge_var]['variable']}").col(f"{cfg_file['variables'][charge_var]['value']}")
        classifier = the_f.get_node(f"/{cfg_file['variables'][classifier_var]['variable']}").col(f"{cfg_file['variables'][classifier_var]['value']}")
        n_gen = cfg_file['juliet'][s][l]['n_files'] * evts_per_file

        weights = weighting.calc_juliet_weight_from_weight_dict_and_prop_matrix(
            weight_dict=weight_dict, prop_matrix=prop_matrix, 
            flux=gzk_partial, n_gen=n_gen, livetime=livetime,
        )

        ehe_weights[s] = np.concatenate((ehe_weights[s], copy.deepcopy(abs(weights))))
        ehe_charge[s] = np.concatenate((ehe_charge[s], copy.deepcopy(charge)))
        ehe_classifier[s] = np.concatenate((ehe_classifier[s], copy.deepcopy(classifier)))

        the_f.close()

ehe_mask = {
    "nue": ehe_charge['nue']>qcut,
    "mu": ehe_charge['mu']>qcut
}


#############################
# muon bundles (corsika)
#############################

cor_classifier = np.asarray([])
cor_charge = np.asarray([])
cor_weights = np.asarray([])

corsika_sets = ["20787"]
for c in corsika_sets:
    print("Working on corsika {}".format(c))
    cor_file = pd.HDFStore(cfg_file['corsika'][c]['file'])
    cor_weighter = weighting.get_weighter( cor_file, 'corsika',
        cfg_file['corsika'][c]['n_files']
        )
    this_cor_charge = cor_weighter.get_column(cfg_file['variables'][charge_var]['variable'], cfg_file['variables'][charge_var]['value'])
    this_cor_classifier = cor_weighter.get_column(cfg_file['variables'][classifier_var]['variable'], cfg_file['variables'][classifier_var]['value'])
    
    L2_q_mask = this_cor_charge > qcut

    cor_classifier = np.concatenate((cor_classifier, this_cor_classifier[L2_q_mask]))
    cor_weights = np.concatenate((cor_weights, cor_weighter.get_weights(cr_flux)[L2_q_mask] * livetime))
    cor_file.close()

if which_method is 'new':
    cor_track_mask = cor_classifier >= classifier_cut
elif which_method is 'old':
    cor_track_mask = cor_classifier <= classifier_cut

cor_weight_total = np.sum(cor_weights)
cor_weight_tracks = np.sum(cor_weights[cor_track_mask])
cor_weight_cascades = np.sum(cor_weights[~cor_track_mask])
cor_misclassifier = cor_weight_cascades/cor_weight_total
print(f"Cor Correct: {cor_weight_tracks/cor_weight_total:.3f}")
print(f"Cor False: {cor_misclassifier:.3f}")



make_plots = True
if make_plots:
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15,5))
    lw = 2

    ax1.hist( 
        ehe_classifier['nue'][ehe_mask['nue']], bins=classifier_bins, 
        weights=ehe_weights["nue"][ehe_mask['nue']],
        histtype='step', label=r'GZK $\nu_{e}$', linewidth=lw)
    n, b, p = ax1.hist( 
        ehe_classifier['mu'][ehe_mask['mu']], bins=classifier_bins, 
        weights=ehe_weights["mu"][ehe_mask['mu']],
        histtype='step', label=r'GZK $\mu$', linewidth=lw)
    ax1.hist( cor_classifier, bins=classifier_bins, weights=cor_weights,
            histtype='step', label='Corsika (H4a)', linewidth=lw)
    ax1.set_yscale('log')
    ax1.set_xlabel(cfg_file['variables'][classifier_var]['label'])
    ax1.set_ylabel('Events / {:.2f} days'.format(livetime/(60*60*24)))
    ax1.set_ylim([1E-6, 1E4])
    ax1.legend()

    ax2.hist( 
        ehe_classifier['nue'][ehe_mask['nue']], bins=classifier_bins, 
        weights=ehe_weights["nue"][ehe_mask['nue']],
        histtype='step', linewidth=lw, density=True)
    ax2.hist( 
        ehe_classifier['mu'][ehe_mask['mu']], bins=classifier_bins, 
        weights=ehe_weights["mu"][ehe_mask['mu']],
        histtype='step', linewidth=lw, density=True)
    ax2.hist( 
        cor_classifier, bins=classifier_bins, 
        weights=cor_weights,
        histtype='step', linewidth=lw, density=True)

    ax2.set_xlabel(cfg_file['variables'][classifier_var]['label'])
    ax2.set_ylabel('Density')

    ax3.hist( 
        ehe_classifier['nue'][ehe_mask['nue']], bins=classifier_bins, 
        weights=ehe_weights["nue"][ehe_mask['nue']],
        histtype='step', linewidth=lw, density=True, cumulative=nue_cumulative_sign)
    ax3.hist( 
        ehe_classifier['mu'][ehe_mask['mu']], bins=classifier_bins, 
        weights=ehe_weights["mu"][ehe_mask['mu']],
        histtype='step', linewidth=lw, density=True, cumulative=-nue_cumulative_sign)
    ax3.hist( 
        cor_classifier, bins=classifier_bins, 
        weights=cor_weights,
        histtype='step', linewidth=lw, density=True, cumulative=-nue_cumulative_sign)
    ax3.set_xlabel(cfg_file['variables'][classifier_var]['label'])
    ax3.set_ylabel('CDF')

    expected_rate = np.sum(ehe_weights['mu'][ehe_mask['mu']])
    print(f"Expected mu rate: {expected_rate}")

    def scaling(ax):
        ax.set_xlim(classifier_lims)
        ax.axvline(classifier_cut, linestyle='--')
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
        scores['nue'][ehe_mask['nue']], 
        scores['mu'][ehe_mask['mu']]
        ))
    y_true = np.concatenate((
        np.full_like(ehe_classifier['nue'], 0)[ehe_mask['nue']],
        np.full_like(ehe_classifier['mu'],  1)[ehe_mask['mu']]
    ))
    y_weights = np.concatenate((
        ehe_weights['nue'][ehe_mask['nue']],
        ehe_weights['mu'][ehe_mask['mu']]
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
