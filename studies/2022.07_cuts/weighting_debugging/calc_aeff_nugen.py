import numpy as np
import tables, yaml, copy
from eheanalysis import weighting, plotting, analysis_9yr
import matplotlib.pyplot as plt
from matplotlib import style
import simweights

style.use('/home/brian/IceCube/ehe/max_tools/EHE_analysis/eheanalysis/ehe.mplstyle')

cfg_file = '../config.yaml'
cfg_file = yaml.safe_load(open(cfg_file))

which_cx = 'cteq5'
qmin = 1E3

neutrino_species = ["nue", "numu", "nutau"]

nugen_aeff = {
    
}

energy_bins, energy_bin_centers = weighting.get_juliet_enu_binning()

for iS, s in enumerate(neutrino_species):

    print("Working on nugen {}".format(s))
    the_f = tables.open_file(cfg_file['nugen'][s]['file'])
    w = simweights.NuGenWeighter(the_f, nfiles = cfg_file['nugen'][s]['n_files'])

    npe = the_f.get_node('/EHEPortiaEventSummarySRT').col('bestNPEbtw')
    fitqual = the_f.get_node('/EHEOpheliaSRT_ImpLF').col('fitQuality')
    recozen = the_f.get_node('/EHEOpheliaParticleSRT_ImpLF').col('zenith')
    
    # muon_bundle_pass = analysis_9yr.muon_bundle_cut_pass_9yr(recozen, npe)
    # track_quality_pass = analysis_9yr.track_quality_cut_pass_9yr(fitqual, npe)
    # total_pass = np.logical_and(muon_bundle_pass, track_quality_pass)
    muon_bundle_pass_3yr = analysis_9yr.muon_bundle_cut_pass_3yr(recozen, npe)
    total_pass = muon_bundle_pass_3yr
      
    ea_avg = w.effective_area(energy_bins, [-1, 1], total_pass)
    
    nugen_aeff[s] = copy.deepcopy(ea_avg)
    
    the_f.close()

import pickle as pickle
output = open('nugen_aeffs.pkl', 'wb')
pickle.dump(
    {
        'energy_bins': energy_bins, 
        'energy_bin_centers': energy_bin_centers,
        'nugen_nu_aeffs': nugen_aeff, 
    },
    output
)

