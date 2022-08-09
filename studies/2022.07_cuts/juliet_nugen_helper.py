import numpy as np
import matplotlib.pyplot as plt
from eheanalysis import weighting, plotting
import copy

def get_juliet_ebin_mask():
    juliet_e_bins, juliet_e_bins_centers = weighting.get_juliet_enu_binning()
    low_mask = np.log10(juliet_e_bins_centers) < 5.5
    high_mask = np.log10(juliet_e_bins_centers) > 7.5
    juliet_surface_mask = np.logical_or(low_mask, high_mask)
    return juliet_surface_mask
    
def mask_and_collapse(to_mask):
    juliet_surface_mask = get_juliet_ebin_mask()
    to_mask_copy = copy.deepcopy(to_mask)
    to_mask_copy[:, juliet_surface_mask] = 0.
    to_mask_copy = to_mask_copy.sum(axis=1) # project down onto the relevant axis
    return to_mask_copy

def make_1D_juliet_nugen_comparison( 
    bins=None, juliet_weights=None, 
    nugen_data=None, nugen_weights=None ):
    
    bin_centers = plotting.get_bin_centers(bins)
    
    fig, (ax, axr) =  plt.subplots(2,1, sharex=True, gridspec_kw={'height_ratios': [3,1]})
    n_j, b_j, p_j = ax.hist(bin_centers, bins=bins, weights=juliet_weights,
        histtype='step', label="Juliet", linewidth=3)
    n_n, b_n, p_n = ax.hist(nugen_data, bins=bins, weights=nugen_weights,
        histtype='step', label="NuGen", linewidth=3, linestyle="--")

    plt.setp(ax.get_xticklabels(), visible=False)
    axr.plot(bin_centers, n_j/n_n - 1, 'o')
    axr.axhline(0, linestyle='--')

    return fig, ax, axr

