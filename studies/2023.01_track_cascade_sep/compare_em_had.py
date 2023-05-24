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
from eheanalysis import weighting, plotting, cuts

style.use('/home/brian/IceCube/ehe/ehe_software/ehe_code/EHE_analysis/eheanalysis/ehe.mplstyle')

from functools import partial
from ehefluxes import fluxes
from eheanalysis import weighting, plotting, cuts
gzk_flux = fluxes.EHEFlux("cosmogenic_ahlers2010_1E18")
gzk_partial = partial(gzk_flux, which_species="nue_sum") 

def astro_flux(energy):
    # flux of mu @ 100 TeV (per-flavor, per-particle flux)
    # then multiply by two to get the nu+nubar sum
    return 2 * 1.44e-18 / 2 * (energy/1e5)**-2.37

# gzk_partial = astro_flux

livetime = 365 * 24 * 60 * 60

cfg_file = 'config.yaml'
cfg_file = yaml.safe_load(open(cfg_file))

which_sample = 'nue_high_energy_with_corrections'

dep_e ={
    "nue": np.asarray([]),
    "mu": np.asarray([])
}
emequiv_e = {
    "nue": np.asarray([]),
    "mu": np.asarray([])
}
int_types = {
    "nue": np.asarray([]),
    "mu": np.asarray([])
}
ehe_weights ={
    "nue": np.asarray([]),
    "mu": np.asarray([])
}
ehe_classifier = {
    "nue": np.asarray([]),
    "mu": np.asarray([])
}

juliet_species = ["nue"]
juliet_energy_levels = ["high_energy"]
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

        the_f = tables.open_file(cfg_file['juliet'][s][l][which_sample]['file'])

        the_depe = the_f.get_node("/DepE").col("value")
        the_emequive = the_f.get_node("/EMEquivVisDepE").col("value")
        the_inttype = the_f.get_node("/NuEProperties").col("GR")
        classifier_var = 'speed'
        classifier = the_f.get_node(f"/{cfg_file['variables'][classifier_var]['variable']}").col(f"{cfg_file['variables'][classifier_var]['value']}")

        weight_dict, prop_matrix, evts_per_file = weighting.get_juliet_weightdict_and_propmatrix(the_f)
        evts_per_file = events_per_file[f"{s}_{l}"] # override to fix L2 issue

        n_gen = 250 * evts_per_file

        weights = weighting.calc_juliet_weight_from_weight_dict_and_prop_matrix(
            weight_dict=weight_dict, prop_matrix=prop_matrix, 
            flux=gzk_partial, n_gen=n_gen, livetime=livetime,
        )
        # weights = np.ones_like(the_inttype)

        dep_e[s] = np.concatenate((dep_e[s], copy.deepcopy(abs(the_depe))))
        emequiv_e[s] = np.concatenate((emequiv_e[s], copy.deepcopy(abs(the_emequive))))
        int_types[s] = np.concatenate((int_types[s], copy.deepcopy(abs(the_inttype))))
        ehe_weights[s] = np.concatenate((ehe_weights[s], copy.deepcopy(abs(weights))))
        ehe_classifier[s] = np.concatenate((ehe_classifier[s], copy.deepcopy(classifier)))
    
        the_f.close()


cmap=plt.cm.plasma

make_plots = True
if make_plots:

    clims = [1E-5, 1E-2]
    norm = colors.LogNorm()
    e_bins = np.logspace(np.log10(1E4),np.log10(1E10),40)
    
    # cut1 = dep_e['nue'] < 1E6
    # cut2 = dep_e['nue'] > 9E5
    # cut = np.logical_and(cut1, cut2)
    # cut = int_types['nue'] > 0
    # cut = dep_e['nue'] >= emequiv_e['nue']
    track_mask = ehe_classifier['nue'] < 0.27
    # cut = np.logical_and(cut, track_mask)
    cut = track_mask

        
    # for i in range(100):
    #     dee = dep_e['nue'][cut][i]
    #     eme = emequiv_e['nue'][cut][i]
    #     print(f"Dep {dee:e}, EM {eme:e}, EM/Dep {eme/dee:.2f}")

    # cut = dep_e['nue'] > 0

    do_e = True
    if do_e:
        
        fig = plt.figure(figsize=(7,5))
        ax = fig.add_subplot(111)
        vals, xedges, yedges, im = plotting.make_2D_hist(
            ax, 
            xvals=dep_e['nue'][cut],
            yvals=emequiv_e['nue'][cut],
            bins=[e_bins, e_bins], cmap=plt.cm.viridis,
            norm=norm,
            weights=ehe_weights['nue'][cut],
            xlabel='Dep E',
            ylabel='EM Equiv E',
        )
        cbar = plt.colorbar(im, ax=ax, label='Evts')
        ax.set_xscale('log')
        ax.set_yscale('log')

        im.set_clim(clims)
        fig.tight_layout()
        ax.set_aspect('equal')
        fig.savefig(f'./figs/em_vs_had.png')
        del fig, ax, im

        # bias
        fig = plt.figure(figsize=(7,5))
        ax = fig.add_subplot(111)
        x, y_med, y_lo, y_hi = plotting.find_contours_2D(
            x_values=dep_e['nue'][cut],
            y_values=(dep_e['nue'][cut] - emequiv_e['nue'][cut])/dep_e['nue'][cut],
            xbins=xedges,
            weights=ehe_weights['nue'][cut],
        )
        ax.plot(x, y_med, color='C0')
        # ax.fill_between(x, y_med, y_hi, color='C0', alpha=0.2)
        # ax.fill_between(x,y_med, y_lo, color='C0', alpha=0.2)
        ax.set_ylabel('(Dep-EM)/Dep')
        ax.set_xlabel('Dep E')
        ax.set_xscale('log')
        # ax.set_ylim([-0.05,0.15])
        ax.set_ylim([-2,2])
        ax.grid()
        ax.hlines(0.,1E6, 1E10, linestyles='--')
        ax.set_xlim(1E4,1E10)
        fig.tight_layout()
        fig.savefig(f'./figs/em_vs_had_bias.png', dpi=300)
        del fig, ax