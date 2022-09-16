import matplotlib
import numpy as np
import tables, yaml, copy
from eheanalysis import weighting, plotting
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import style
import pandas as pd
style.use('/home/brian/IceCube/ehe/max_tools/EHE_analysis/eheanalysis/ehe.mplstyle')

cfg_file = '../config.yaml'
cfg_file = yaml.safe_load(open(cfg_file))

juliet_species = ["nue", "numu", "nutau", "mu", "tau"]
juliet_energy_levels = ["high_energy", "very_high_energy"]
# juliet_species = ["tau"]
# juliet_energy_levels = ["very_high_energy"]

fig, ax = plt.subplots(1, 1, figsize=(6,5))
bins = np.logspace(3, 12, 31)

for j in juliet_species:
    ehe_energies = np.asarray([])
    for l in juliet_energy_levels:
        print(f"Working on {j}, {l}")
        the_file = cfg_file['juliet'][j][l]['file']
        with tables.open_file(the_file) as f:
            juliet_primary = f.get_node('/I3JulietPrimaryParticle').col('energy')
            ehe_energies = np.concatenate((ehe_energies, juliet_primary))
    ax.hist(ehe_energies, bins=bins, 
            histtype='step', linewidth=3, label="{}".format(j))
        
ax.set_xlabel("Primary Energy / Gev")
ax.set_ylabel("Unweighted Counts")
ax.set_xscale('log')
ax.set_yscale('log')
ax.legend(loc='upper left')
fig.tight_layout()
fig.savefig("unweighted_juliet_counts.png")