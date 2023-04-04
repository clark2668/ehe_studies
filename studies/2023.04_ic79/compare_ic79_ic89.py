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
from eheanalysis import plotting, cuts


livetime = 365*24*60*60

cfg_file = 'config.yaml'
cfg_file = yaml.safe_load(open(cfg_file))

def harvest_values(cfg_file, dataset):
    return_dict = {}
    the_f = tables.open_file(cfg_file['burn_sample'][dataset]['file'])
    charge = the_f.get_node(f"/{cfg_file['variables']['charge']['variable']}").col(f"{cfg_file['variables']['charge']['value']}")
    zen = the_f.get_node(f"/{cfg_file['variables']['zenith']['variable']}").col(f"{cfg_file['variables']['zenith']['value']}")
    speed = the_f.get_node(f"/{cfg_file['variables']['speed']['variable']}").col(f"{cfg_file['variables']['speed']['value']}")
    ndoms = the_f.get_node(f"/{cfg_file['variables']['ndoms']['variable']}").col(f"{cfg_file['variables']['ndoms']['value']}")
    return_dict['charge'] = copy.deepcopy(charge)
    return_dict['zenith'] = copy.deepcopy(zen)
    return_dict['speed'] = copy.deepcopy(speed)
    return_dict['ndoms'] = copy.deepcopy(ndoms)
    return_dict['livetime'] = cfg_file['burn_sample'][dataset]['livetime']
    the_f.close()
    
    L1_mask = return_dict['charge'] > 1E3

    q_mask = return_dict['charge'] > 27500
    ndoms_mask = return_dict['ndoms'] > 100
    L2_mask = np.logical_and(q_mask, ndoms_mask)

    speed_mask = cuts.track_quality_cut(return_dict['speed'], return_dict['charge'])
    L3_mask = np.logical_and(L2_mask, speed_mask)
    
    return_dict['L1_mask'] = L1_mask
    return_dict['L2_mask'] = L2_mask
    return_dict['L3_mask'] = L3_mask
    return return_dict

ic79_dict = harvest_values(cfg_file, 'IC79-2010-pass2')
ic86_dict = harvest_values(cfg_file, 'IC86-2011-pass2')

def make_compare_plot(var_bins, ic79_vals, ic86_vals):

    fig, (ax, axr) = plt.subplots(2, 1, figsize=(5,5), sharex=True, gridspec_kw={'height_ratios': [3, 1]})
    bin_centers = plotting.get_bin_centers(var_bins)

    ic79_sum, b = np.histogram(ic79_vals,bins=var_bins)
    ic86_sum, b = np.histogram(ic86_vals,bins=var_bins)

    ic79_errs = np.sqrt(ic79_sum)
    ic86_errs = np.sqrt(ic86_sum)

    ic79_rate = ic79_sum/ic79_dict['livetime']
    ic86_rate = ic86_sum/ic86_dict['livetime']

    ic79_rate_errs = ic79_errs/ic79_dict['livetime']
    ic86_rate_errs = ic86_errs/ic86_dict['livetime']

    do_rate = True
    if do_rate:
        which_to_plot_ic79 = ic79_rate
        which_to_plot_ic86 = ic86_rate
        which_errs_to_plot_ic79 = ic79_rate_errs
        which_errs_to_plot_ic86 = ic86_rate_errs
    else:
        which_to_plot_ic79 = ic79_sum
        which_to_plot_ic86 = ic86_sum
        which_errs_to_plot_ic79 = ic79_errs
        which_errs_to_plot_ic86 = ic86_errs

    do_errors = True
    if do_errors:

        ax.errorbar(
            x=bin_centers,
            y=which_to_plot_ic79,
            yerr=which_errs_to_plot_ic79,
            fmt='o', alpha=0.5,
            label=f'IC79, {np.sum(which_to_plot_ic79):.2e} Hz'
        )
        ax.errorbar(
            x=bin_centers,
            y=which_to_plot_ic86,
            yerr=which_errs_to_plot_ic86,
            fmt='s', alpha=0.5,
            label=f'IC86, {np.sum(which_to_plot_ic86):.2e} Hz'
        )

        ratio = which_to_plot_ic79/which_to_plot_ic86

        err_term_a = (which_errs_to_plot_ic79**2)/(which_to_plot_ic79**2)
        err_term_b = (which_errs_to_plot_ic86**2)/(which_to_plot_ic86**2)
        err_term_tot = np.sqrt(err_term_a + err_term_b)
        err_term_tot = ratio * err_term_tot

        overall_rate_diff = np.sum(which_to_plot_ic79)/np.sum(which_to_plot_ic86)

        axr.errorbar(
            x=bin_centers,
            y=ratio,
            yerr=err_term_tot,
            fmt='ko'
        )

    else:
        ic79_sum, b, p = ax.hist(
            x=bin_centers, weights=which_to_plot_ic79, bins=var_bins,
            label=f'IC79, {np.sum(which_to_plot_ic79):.2e} Hz', histtype='step', lw=2
        )

        ic86_sum, b, p = ax.hist(
            x=bin_centers, weights=which_to_plot_ic86, bins=var_bins, 
            label=f'IC86, {np.sum(which_to_plot_ic86):.2e} Hz', histtype='step', lw=2
        )

        axr.plot(bin_centers, ic79_sum/ic86_sum,'o-')

    ax.legend()
    ax.set_ylabel("Hz")
    axr.set_ylabel('IC79/IC86')
    axr.set_ylim([0.5,1.5])
    axr.grid()
    axr.axhline(y=overall_rate_diff,linestyle='--', color='red',
                label=f'Rate Ratio: {overall_rate_diff:.2f}'
                )
    axr.legend()
    return fig, ax, axr

levels = ['L1', 'L2', 'L3']
for level in levels:

    def set_yscale(ax, axr, level):

        if level is 'L1':
            ax.set_ylim([1E-7, 1E-2])
        if level is 'L2':
            ax.set_ylim([1E-7, 1E-4])
        if level is 'L3':
            ax.set_ylim([1E-7, 2E-5])


    q_bins = np.logspace(4,6,20)
    fig, ax, axr = make_compare_plot(q_bins, 
                                    ic79_dict['charge'][ic79_dict[f'{level}_mask']],
                                    ic86_dict['charge'][ic86_dict[f'{level}_mask']]
                                    )
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_title(f"After {level}")
    set_yscale(ax, axr, level)
    axr.set_xlabel(r'HQTot')
    fig.tight_layout()
    fig.savefig(f'plots/charge_rate_{level}.png', dpi=300)
    del fig, ax, axr

    czen_bins = np.linspace(-1,1,20)
    fig, ax, axr = make_compare_plot(czen_bins, 
                                    np.cos(ic79_dict['zenith'][ic79_dict[f'{level}_mask']]),
                                    np.cos(ic86_dict['zenith'][ic86_dict[f'{level}_mask']])
                                    )
    ax.set_yscale('log')
    set_yscale(ax, axr, level)
    ax.set_title(f"After {level}")
    axr.set_xlabel(r'cos($\theta$)')
    fig.tight_layout()
    fig.savefig(f'plots/czen_rate_{level}.png', dpi=300)
    del fig, ax, axr