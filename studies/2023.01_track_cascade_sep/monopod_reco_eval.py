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
from eheanalysis import weighting, plotting, cuts
gzk_flux = fluxes.EHEFlux("cosmogenic_ahlers2010_1E18")
gzk_partial = partial(gzk_flux, which_species="nue_sum") 

style.use('/home/brian/IceCube/ehe/ehe_software/ehe_code/EHE_analysis/eheanalysis/ehe.mplstyle')

livetime = 365 * 24 * 60 * 60
print(livetime)

cfg_file = 'config.yaml'
cfg_file = yaml.safe_load(open(cfg_file))

qcut = 27500
charge_var = 'hqtot'
classifier_var = 'speed'
classifier_bins = np.linspace(0, 0.5, 150)
classifier_cut = 0.27
ndoms_cut = 100

ehe_weights ={
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
ehe_splinempe_zen = {
    "nue": np.asarray([]),
    "mu": np.asarray([])
}
ehe_splinempe_azi = {
    "nue": np.asarray([]),
    "mu": np.asarray([])
}
ehe_monopod_zen = {
    "nue": np.asarray([]),
    "mu": np.asarray([])
}
ehe_monopod_azi = {
    "nue": np.asarray([]),
    "mu": np.asarray([])
}
ehe_true_zen = {
    "nue": np.asarray([]),
    "mu": np.asarray([])
}
ehe_true_azi = {
    "nue": np.asarray([]),
    "mu": np.asarray([])
}



juliet_species = ["nue"]
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
        ndoms = the_f.get_node(f"/{cfg_file['variables']['ndoms']['variable']}").col(f"{cfg_file['variables']['ndoms']['value']}")
        classifier = the_f.get_node(f"/{cfg_file['variables'][classifier_var]['variable']}").col(f"{cfg_file['variables'][classifier_var]['value']}")
        splinempe_zen = the_f.get_node("/EHE_SplineMPE").col("zenith")
        splinempe_azi = the_f.get_node("/EHE_SplineMPE").col("azimuth")
        monopod_zen = the_f.get_node("/EHE_Monopod").col("zenith")
        monopod_azi = the_f.get_node("/EHE_Monopod").col("azimuth")
        truth_zen = the_f.get_node("/PolyplopiaPrimary").col("zenith")
        truth_azi = the_f.get_node("/PolyplopiaPrimary").col("azimuth")

        n_gen = cfg_file['juliet'][s][l]['n_files'] * evts_per_file

        weights = weighting.calc_juliet_weight_from_weight_dict_and_prop_matrix(
            weight_dict=weight_dict, prop_matrix=prop_matrix, 
            flux=gzk_partial, n_gen=n_gen, livetime=livetime,
        )

        ehe_weights[s] = np.concatenate((ehe_weights[s], copy.deepcopy(abs(weights))))
        ehe_charge[s] = np.concatenate((ehe_charge[s], copy.deepcopy(charge)))
        ehe_ndoms[s] = np.concatenate((ehe_ndoms[s], copy.deepcopy(ndoms)))
        ehe_classifier[s] = np.concatenate((ehe_classifier[s], copy.deepcopy(classifier)))
        
        ehe_splinempe_zen[s] = np.concatenate((ehe_splinempe_zen[s], copy.deepcopy(splinempe_zen)))
        ehe_splinempe_azi[s] = np.concatenate((ehe_splinempe_azi[s], copy.deepcopy(splinempe_azi)))
        
        ehe_monopod_zen[s] = np.concatenate((ehe_monopod_zen[s], copy.deepcopy(monopod_zen)))
        ehe_monopod_azi[s] = np.concatenate((ehe_monopod_azi[s], copy.deepcopy(monopod_azi)))
        
        ehe_true_zen[s] = np.concatenate((ehe_true_zen[s], copy.deepcopy(truth_zen)))
        ehe_true_azi[s] = np.concatenate((ehe_true_azi[s], copy.deepcopy(truth_azi)))

        the_f.close()

ehe_mask = {
    "nue": ehe_charge['nue']>qcut,
    "mu": ehe_charge['mu']>qcut
}

# build masks
for f in ehe_mask.keys():
    q_mask = ehe_charge[f] > qcut
    ndoms_mask = ehe_ndoms[f] > ndoms_cut
    track_qual_mask = cuts.track_quality_cut(ehe_classifier[f], ehe_charge[f])
    track_mask = ehe_classifier[f] < 0.27
    total_mask = np.logical_and(q_mask, ndoms_mask)
    total_mask = np.logical_and(total_mask, track_qual_mask)
    total_mask = np.logical_and(total_mask, track_mask)
    ehe_mask[f] = total_mask
    # ehe_mask[f] = ehe_charge[f] > 0

cmap=plt.cm.plasma

make_plots = True
if make_plots:
    fig, axarr  = plt.subplots(2, 3, figsize=(15,7))
    axarr = axarr.flatten()

    czen_bins = np.linspace(-1, 1, 40)
    azi_bins = np.linspace(0,np.pi*2, 40)
    clims = [1E-5, 1E-2]

    a, b, c, im0 = axarr[0].hist2d(
        np.cos(ehe_true_zen['nue'][ehe_mask['nue']]), np.cos(ehe_splinempe_zen['nue'][ehe_mask['nue']]),
        bins=[czen_bins, czen_bins],
        weights=ehe_weights['nue'][ehe_mask['nue']],
        norm=colors.LogNorm(),
    )
    axarr[0].set_xlabel('True czen'); axarr[0].set_ylabel('SplineMPE czen')
    cbar0 = plt.colorbar(im0, ax=axarr[0], label='Evts/Year')
    im0.set_clim(clims)

    a, b, c, im1 = axarr[1].hist2d(
        np.cos(ehe_true_zen['nue'][ehe_mask['nue']]), np.cos(ehe_monopod_zen['nue'][ehe_mask['nue']]),
        bins=[czen_bins, czen_bins],
        weights=ehe_weights['nue'][ehe_mask['nue']],
        norm=colors.LogNorm(),
    )
    axarr[1].set_xlabel('True czen'); axarr[1].set_ylabel('Monopod czen')
    cbar1 = plt.colorbar(im1, ax=axarr[1], label='Evts/Year')
    im1.set_clim(clims)

    a, b, c, im3 = axarr[3].hist2d(
        ehe_true_azi['nue'][ehe_mask['nue']], ehe_splinempe_azi['nue'][ehe_mask['nue']],
        bins=[azi_bins, azi_bins],
        weights=ehe_weights['nue'][ehe_mask['nue']],
        norm=colors.LogNorm(),
    )
    axarr[3].set_xlabel('True azi'); axarr[3].set_ylabel('SplineMPE azi')
    cbar3 = plt.colorbar(im1, ax=axarr[3], label='Evts/Year')
    im3.set_clim(clims)

    a, b, c, im4 = axarr[4].hist2d(
        ehe_true_azi['nue'][ehe_mask['nue']], ehe_monopod_azi['nue'][ehe_mask['nue']],
        bins=[azi_bins, azi_bins],
        weights=ehe_weights['nue'][ehe_mask['nue']],
        norm=colors.LogNorm(),
    )
    axarr[4].set_xlabel('True azi'); axarr[4].set_ylabel('Monopod azi')
    cbar4 = plt.colorbar(im1, ax=axarr[4], label='Evts/Year')
    im4.set_clim(clims)

    fig.tight_layout()
    fig.savefig(f'./figs/reco_comparison.png')


