import numpy as np
import pickle
from copy import copy
from skimage.measure import block_reduce
from sklearn.metrics import roc_curve
from sklearn.metrics import auc
from scipy.optimize import curve_fit
from scipy.optimize import brentq
from scipy.signal import argrelmin
from scipy.integrate import quad

import sys
import os
if os.uname().nodename in ['Brians-MBP', 'dhcp-172-16-174-138.private.pa.umd.edu']:
    sys.path.append('/Users/brianclark/Documents/work/IceCube/ehe/ehe_deps/energy_loss_pdfs')
elif os.uname().nodename in ['condor00']:
    sys.path.append('/data/i3home/baclark/IceCube/ehe/career_plots/tools/energy_loss_pdfs')
else:
    raise NotImplementedError

from likelihood import Likelihood_1D


def exp(x, a, b):
    return a * np.exp(b * x)


def zombie_pdf(spline, exp_func):
    elost_vals = np.logspace(-10, 0, 1001)
    spline_max = np.argmax(spline(elost_vals))

    def shifted_spline(x, shift=1e-5):
        return spline(x) - shift

    rel_mins = argrelmin(spline(elost_vals[:spline_max]), order=3)
    left_pivot = brentq(shifted_spline, elost_vals[np.max(rel_mins)], elost_vals[spline_max])
    left_idx = np.argmin(np.abs(elost_vals - left_pivot))
    popt, pcov = curve_fit(
        exp_func,
        elost_vals[left_idx:left_idx+10],
        spline(elost_vals[left_idx:left_idx+10]),
        maxfev=10000)

    rel_mins = argrelmin(spline(elost_vals[spline_max:]), order=3) + spline_max
    right_pivot = brentq(shifted_spline, elost_vals[spline_max], elost_vals[rel_mins[0][1]])

    right_idx = np.argmin(np.abs(elost_vals - right_pivot))
    popt_, pcov_ = curve_fit(
        exp_func,
        elost_vals[right_idx-20:right_idx],
        spline(elost_vals[right_idx-20:right_idx]),
        maxfev=10000)

    def func(x):
        if not isinstance(x, np.ndarray):
            x = np.array([x])
        values = np.zeros_like(x)
        mask = x < left_pivot * 1.15
        values[mask] = exp_func(x[mask], *popt)
        values[~mask] = spline(x[~mask])
        mask = x > right_pivot * 0.85
        values[mask] = 1e-10
        values[mask] = exp_func(x[mask], *popt_)
        return values

    return func


def get_1d_rlogl(
        mc_df,
        likelihood,
        table_name='EHEMuMillipede_SplineMPEseed_1e-12_vecd',
        key_name='rlogl_1d_40m_avg',
        min_bins=30,
        n_bins_to_combine=4,
        use_average_instead_of_median=False,
        return_auc=False):

    # Get all relevant columns from the dataframe
    e_deposition_cols = [
        col for col in mc_df.columns
        if col.startswith(table_name + '.vector_elem_')]
    e_depositions = np.array(mc_df[e_deposition_cols])

    # Combine n adjacent bins to reduce reconstruction smearing
    e_depositions = block_reduce(e_depositions,
                                 block_size=(1, n_bins_to_combine),
                                 func=np.nanmean) * n_bins_to_combine
    # Mask out skimming events with less than min_bins inside
    # the detector volume or only (0, nan) energy depositions
    non_skimming = np.sum(np.isfinite(e_depositions), axis=1) >= min_bins
    non_skimming = np.logical_and(non_skimming,
                                  np.nansum(e_depositions, axis=1) > 0)
    # mc_df['not_skimming'] = False
    # mc_df.loc[non_skimming, 'not_skimming'] = True

    e_depositions_ = copy(e_depositions)
    elost_center = np.logspace(-10, 0, 1001)

    if not use_average_instead_of_median:
        spline_cumsum = np.cumsum(likelihood.pdf(elost_center))
        idx = np.searchsorted(spline_cumsum, spline_cumsum[-1] * 0.5)
        pdf_median = elost_center[idx]
        reco_median = np.nanmedian(e_depositions_, axis=1)

        # Scale the reconstructed energy losses and calculate the rlogl values
        rel_e_depositions = e_depositions_ / reco_median[:, np.newaxis] * pdf_median
    else:
        reco_average = np.nanmean(e_depositions_, axis=1)
        pdf_norm, _ = quad(lambda x: spline(x), 1e-10, 1)
        pdf_average, _ = quad(lambda x: x * spline(x) / pdf_norm, 1e-10, 1)

        rel_e_depositions = e_depositions_ / reco_average[:, np.newaxis] * pdf_average

    rlogl = likelihood.rlogl(rel_e_depositions[non_skimming])

    # Check the fraction of llh values that are nan/inf
    finite = np.isfinite(rlogl)
    # print(np.sum(finite), np.sum(~finite))

    # sig = mc_df.loc[non_skimming, 'ROC_Label'] == 1
    # finite_sig = np.logical_and(finite, sig)
    # not_finite_sig = np.logical_and(~finite, sig)
    # bkg = mc_df.loc[non_skimming, 'ROC_Label'] == 0
    # finite_bkg = np.logical_and(finite, bkg)
    # not_finite_bkg = np.logical_and(~finite, bkg)
    # print(np.sum(sig), np.sum(not_finite_sig))
    # print(np.sum(mc_df.loc[non_skimming, 'Weight'][sig]), np.sum(mc_df.loc[non_skimming, 'Weight'][not_finite_sig]))
    # print(np.sum(mc_df.loc[non_skimming, 'Weight'][bkg]), np.sum(mc_df.loc[non_skimming, 'Weight'][not_finite_bkg]))
    # print(np.sum(bkg), np.sum(not_finite_bkg))

    mc_df[key_name] = np.nan
    mc_df.loc[non_skimming, key_name] = rlogl

    if return_auc:
        from neutrino_level.steps.general_modules.label_maker import add_label_from_components_to_df
        add_label_from_components_to_df(
            mc_df,
            colname='ROC_Label',
            sig_comps=[2],
            bkg_comps=[6]
        )
        mask = np.isfinite(mc_df[key_name])
        fpr, tpr, threshs = roc_curve(
            mc_df.loc[mask, 'ROC_Label'],
            mc_df.loc[mask, key_name],
            sample_weight=mc_df.loc[mask, 'Weight']
        )
        return auc(fpr, tpr)
    else:
        return None


def get_1d_rlogl_inject_mip_losses(
        mc_df,
        likelihood,
        table_name='EHEMuMillipede_SplineMPEseed_1e-12_vecd',
        key_name='rlogl_1d_40m_avg',
        min_bins=30,
        n_bins_to_combine=4,
        return_auc=False):

    # Get all relevant columns from the dataframe
    e_deposition_cols = [
        col for col in mc_df.columns
        if col.startswith(table_name + '.vector_elem_')]
    e_depositions = np.array(mc_df[e_deposition_cols])

    # Combine n adjacent bins to reduce reconstruction smearing
    e_depositions = block_reduce(e_depositions,
                                 block_size=(1, n_bins_to_combine),
                                 func=np.nanmean) * n_bins_to_combine
    # Mask out skimming events with less than min_bins inside
    # the detector volume or only (0, nan) energy depositions
    non_skimming = np.sum(np.isfinite(e_depositions), axis=1) >= min_bins
    non_skimming = np.logical_and(non_skimming,
                                  np.nansum(e_depositions, axis=1) > 0)
    # mc_df['not_skimming'] = False
    # mc_df.loc[non_skimming, 'not_skimming'] = True

    e_depositions_ = copy(e_depositions)

    ionization_loss = 2.32467146 * n_bins_to_combine
    # Alternatively replace all zeros with ionization loss
    if ionization_loss is None:
        raise ValueError(
            'ionization_loss cant be None '
            'when zeroes are not ignored.'
        )
    msk = np.cumsum(e_depositions_, axis=1) == 0
    e_depositions_[msk] = np.nan
    e_depositions_[e_depositions_ == 0] = ionization_loss

    elost_center = np.logspace(-10, 0, 1001)
    spline_cumsum = np.cumsum(likelihood.pdf(elost_center))
    idx = np.searchsorted(spline_cumsum, spline_cumsum[-1] * 0.5)
    pdf_median = elost_center[idx]
    reco_median = np.nanmedian(e_depositions_, axis=1)

    # Scale the reconstructed energy losses and calculate the rlogl values
    rel_e_depositions = e_depositions_ / reco_median[:, np.newaxis] * pdf_median
    rlogl = likelihood.rlogl(rel_e_depositions[non_skimming])

    # Check the fraction of llh values that are nan/inf
    finite = np.isfinite(rlogl)
    # print(np.sum(finite), np.sum(~finite))

    # sig = mc_df.loc[non_skimming, 'ROC_Label'] == 1
    # finite_sig = np.logical_and(finite, sig)
    # not_finite_sig = np.logical_and(~finite, sig)
    # bkg = mc_df.loc[non_skimming, 'ROC_Label'] == 0
    # finite_bkg = np.logical_and(finite, bkg)
    # not_finite_bkg = np.logical_and(~finite, bkg)
    # print(np.sum(sig), np.sum(not_finite_sig))
    # print(np.sum(mc_df.loc[non_skimming, 'Weight'][sig]), np.sum(mc_df.loc[non_skimming, 'Weight'][not_finite_sig]))
    # print(np.sum(mc_df.loc[non_skimming, 'Weight'][bkg]), np.sum(mc_df.loc[non_skimming, 'Weight'][not_finite_bkg]))
    # print(np.sum(bkg), np.sum(not_finite_bkg))

    mc_df[key_name] = np.nan
    mc_df.loc[non_skimming, key_name] = rlogl

    if return_auc:
        from neutrino_level.steps.general_modules.label_maker import add_label_from_components_to_df
        add_label_from_components_to_df(
            mc_df,
            colname='ROC_Label',
            sig_comps=[2],
            bkg_comps=[6]
        )
        mask = np.isfinite(mc_df[key_name])
        fpr, tpr, threshs = roc_curve(
            mc_df.loc[mask, 'ROC_Label'],
            mc_df.loc[mask, key_name],
            sample_weight=mc_df.loc[mask, 'Weight']
        )
        return auc(fpr, tpr)
    else:
        return None


def get_1d_rlogl_with_deposited_energy(
        mc_df,
        likelihood,
        table_name='EHEMuMillipede_SplineMPEseed_1e-12_vecd',
        key_name='rlogl_1d_40m_avg',
        min_bins=30,
        n_bins_to_combine=4,
        return_auc=False):

    # Get all relevant columns from the dataframe
    e_deposition_cols = [
        col for col in mc_df.columns
        if col.startswith(table_name + '.vector_elem_')]
    e_depositions = np.array(mc_df[e_deposition_cols])

    # Combine n adjacent bins to reduce reconstruction smearing
    e_depositions = block_reduce(e_depositions,
                                 block_size=(1, n_bins_to_combine),
                                 func=np.nanmean) * n_bins_to_combine
    # Mask out skimming events with less than min_bins inside
    # the detector volume or only (0, nan) energy depositions
    non_skimming = np.sum(np.isfinite(e_depositions), axis=1) >= min_bins
    non_skimming = np.logical_and(non_skimming,
                                  np.nansum(e_depositions, axis=1) > 0)

    e_depositions_ = copy(e_depositions)
    reco_deposited = np.nansum(e_depositions_, axis=1)

    # Scale the reconstructed energy losses and calculate the rlogl values
    rel_e_depositions = e_depositions_ / reco_deposited[:, np.newaxis]
    rlogl = likelihood.rlogl(rel_e_depositions[non_skimming])

    # Check the fraction of llh values that are nan/inf
    finite = np.isfinite(rlogl)
    # print(np.sum(finite), np.sum(~finite))

    mc_df[key_name] = np.nan
    mc_df.loc[non_skimming, key_name] = rlogl

    if return_auc:
        from neutrino_level.steps.general_modules.label_maker import add_label_from_components_to_df
        add_label_from_components_to_df(
            mc_df,
            colname='ROC_Label',
            sig_comps=[2],
            bkg_comps=[6]
        )
        mask = np.isfinite(mc_df[key_name])
        fpr, tpr, threshs = roc_curve(
            mc_df.loc[mask, 'ROC_Label'],
            mc_df.loc[mask, key_name],
            sample_weight=mc_df.loc[mask, 'Weight']
        )
        return auc(fpr, tpr)
    else:
        return None


if os.uname().nodename in ['Brians-MBP', 'dhcp-172-16-174-138.private.pa.umd.edu']:
    filename = '/Users/brianclark/Documents/work/IceCube/ehe/ehe_deps/energy_loss_pdfs/pdf_slice_spline_40m.pkl'
elif os.uname().nodename in ['condor00']:
    filename = '/data/i3home/baclark/IceCube/ehe/career_plots/tools/energy_loss_pdfs/pdf_slice_spline_40m.pkl'
else:
    raise NotImplementedError

with open(filename, 'rb') as open_file:
    spline = pickle.load(open_file)

stitched_spline = zombie_pdf(spline, exp)

likelihood = Likelihood_1D(stitched_spline)

