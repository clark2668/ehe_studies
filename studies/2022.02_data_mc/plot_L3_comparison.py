import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import argparse

import utils_weights, utils_ehe
livetime = 86400 * 365
livetime = 350 # in seconds

parser = argparse.ArgumentParser()
parser.add_argument("-numu", type=str,
    dest="numu_file", required=True,
    help="paths to numu file")
parser.add_argument("-nue", type=str, 
    dest="nue_file", required=True,
    help="paths to nue file")
parser.add_argument("-cor", type=str, 
    dest="cor_file", required=True,
    help="paths to corsika file")
parser.add_argument("-data", type=str,
    dest="data_file", required=True,
    help="path to the 10% data file")
args = parser.parse_args()

numu_file = pd.HDFStore(args.numu_file)
nue_file = pd.HDFStore(args.nue_file)
cor_file = pd.HDFStore(args.cor_file)

which_one = 'new'
if which_one is 'original':
    charge_var = ['EHEPortiaEventSummarySRT', 'bestNPEbtw']
    zenith_var = ['EHEOpheliaParticleSRT_ImpLF', 'zenith']
    fit_var = ['EHEOpheliaSRT_ImpLF', 'fitQuality']
    charge_label = "Portia"
    zenith_label = "Ophelia"
elif which_one is 'new':
    charge_var = ['Homogenized_QTot', 'value']
    zenith_var = ['LineFit', 'zenith']
    fit_var = ['LineFitQuality', 'value']
    charge_label = "HQtot"
    zenith_label = "LineFit"

atmo_model = 'H3a_SIBYLL23C'
cr_model = 'GaisserH3a'
atmo_flux_model = utils_weights.get_flux_model(atmo_model, 'nugen')
cr_flux_model = utils_weights.get_flux_model(cr_model, 'corsika')

numu_weighter = utils_weights.get_weighter(numu_file, 'nugen', 1000)
nue_weighter = utils_weights.get_weighter(nue_file, 'nugen', 1000)
cor_weighter = utils_weights.get_weighter(cor_file, 'corsika', 1000)

# # numu_zenith = numu_weighter.get_column(zenith_var[0], zenith_var[1])
# # numu_chisqured = numu_weighter.get_column(fit_var[0], fit_var[1])
numu_npe = numu_weighter.get_column(charge_var[0], charge_var[1])

# # nue_zenith = nue_weighter.get_column(zenith_var[0], zenith_var[1])
# # nue_chisqured = nue_weighter.get_column(fit_var[0], fit_var[1])
nue_npe = nue_weighter.get_column(charge_var[0], charge_var[1])

# # cor_zenith = cor_weighter.get_column(zenith_var[0], zenith_var[1])
# # cor_chisqured = cor_weighter.get_column(fit_var[0], fit_var[1])
cor_npe = cor_weighter.get_column(charge_var[0], charge_var[1])

numu_atmo_weights = numu_weighter.get_weights(atmo_flux_model)
nue_atmo_weights = nue_weighter.get_weights(atmo_flux_model)
cor_weights = cor_weighter.get_weights(cr_flux_model)

numu_file.close()
nue_file.close()
cor_file.close()

numu_atmo_weights *= livetime
nue_atmo_weights *= livetime
cor_weights *= livetime

data_file = pd.HDFStore(args.data_file)
values = np.asarray(data_file.get('Homogenized_QTot').get('value'))
data_npe = np.asarray(data_file.get(charge_var[0]).get(charge_var[1]))
data_npe = np.log10(data_npe)

do_L2_plot = True
if do_L2_plot:

    fig = plt.figure()
    ax = fig.add_subplot(111)
    bins = np.linspace(2, 6, 40)

    # histogram and then plot the data
    data, data_bins = np.histogram(data_npe, bins=bins)
    errs = np.sqrt(data)
    binscenters = np.array([0.5 * (data_bins[i] + data_bins[i+1]) for i in range(len(data_bins)-1)])
    ax.errorbar(binscenters, data, yerr=errs, fmt='ko', label='Burn Sample')

    # histogram the backgrounds
    ax.hist(
        [np.log10(numu_npe), np.log10(nue_npe), np.log10(cor_npe)],
        weights=[numu_atmo_weights, nue_atmo_weights, cor_weights],
        label=[
            r'Atm $\nu_{\mu}$' + ', {}'.format(atmo_model), 
            r'Atm $\nu_{e}$'+ ', {}'.format(atmo_model), 
            r'Atm $\mu$' + ', {}'.format(cr_model)
            ],
        bins=bins,
        stacked=True
    )

    ax.set_ylabel('Events (livetime ~350 seconds)')
    ax.set_xlabel('Charge ({}, NPE)'.format(charge_label))
    ax.set_yscale('log')
    ax.set_ylim([1E-7, 1E3])
    ax.legend()
    plt.tight_layout()
    fig.savefig('burn_sample_demo.png')

data_file.close()