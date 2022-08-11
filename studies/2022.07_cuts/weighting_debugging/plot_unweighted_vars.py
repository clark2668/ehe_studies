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

bins_czen = np.linspace(-1, 1, 20)
bins_azi = np.linspace(0, 2*np.pi, 20)

juliet_species = ["nue", "numu", "nutau", "mu", "tau"]

geo_info = np.load("dom_loc.npz")
om_x = geo_info["om_x"]
om_z = geo_info["om_z"]
strings_x = geo_info["strings_x"]
strings_y = geo_info["strings_y"]


for j in juliet_species:
    print(f"Working on {j}")
    the_file = cfg_file['juliet'][j]['high_energy']['file_cteq5']
    with tables.open_file(the_file) as f:
        juliet_primary = f.get_node('/I3JulietPrimaryParticle')
        first_mctree = f.get_node('/PrimaryEvent')
        impact_param = f.get_node('/ClosestApproach').col('value')
        # x, y = np.unique(first_mctree.col('pdg_encoding'), return_counts=True)
        # print(y)
        
        fig, (ax1, ax2) = plt.subplots(1,2, figsize=(10,5))
        evs1,_ ,__ = ax1.hist(np.cos(juliet_primary.col('zenith')), bins=bins_czen, 
                             histtype='step', linewidth=3, label='I3JulietPrimaryParticle')
        evs2,_ ,__ = ax1.hist(np.cos(first_mctree.col('zenith')), bins=bins_czen, 
                             histtype='step', linewidth=3, label='I3MCTree Primary', linestyle='--')
        ax1.set_xlabel(r'cos($\theta$)')
        ax1.set_ylabel('Counts')
        ax1.legend(loc='lower right')

        evs1,_ ,__ = ax2.hist(juliet_primary.col('azimuth'), bins=bins_azi, 
                             histtype='step', linewidth=3, label='I3JulietPrimaryParticle')
        evs1,_ ,__ = ax2.hist(first_mctree.col('azimuth'), bins=bins_azi, 
                             histtype='step', linewidth=3, label='I3MCTree Primary', linestyle='--')
        ax2.set_xlabel(r'$\phi$')
        ax2.set_ylabel('Counts')
        fig.suptitle(f"{j}")
        plt.tight_layout()
        fig.savefig(f'./unweighted_plots/unweighted_angles_{j}.png')
        del fig, ax1, ax2
        
        fig, ax = plt.subplots(1,1)
        bins = np.linspace(0, 900, 20)
        ax.hist(impact_param, histtype='step', linewidth=3, bins=bins)
        ax.set_xlabel('Distance of Closest Approach to I3 Center / m ')
        ax.set_ylabel("Unweighted Events")
        ax.set_title(f"{j}")
        ax.set_ylim([0, 16000])
        plt.tight_layout()
        fig.savefig(f'./unweighted_plots/unweighted_impactparam_{j}.png')
        del fig, ax
        
        bins = [np.linspace(-1000, 1000, 51), np.linspace(-1000, 1000, 51)]
        
        fig = plt.figure(figsize=(15,5))
        ax1 = fig.add_subplot(131)

        evs1, _, _, im1  = ax1.hist2d(first_mctree.col('x'),first_mctree.col('y'),
                                      bins=bins, cmin=1, norm=colors.LogNorm())
        ax1.set_xlabel("X")
        ax1.set_ylabel("Y")
        ax1.plot(strings_x, strings_y, 'x', color='r', markersize=2)

        c1 = plt.Circle((0, 0), 880, color='r', fill=False, linestyle='--')
        plt.gca().add_patch(c1)
        cbar1 = plt.colorbar(im1, ax=ax1)
        cbar1.set_label('Counts')
        ax1.set_aspect('equal')

        ax2 = fig.add_subplot(132)
        evs2, _, _, im2  = ax2.hist2d(first_mctree.col('x'),first_mctree.col('z'),
                                      bins=bins, cmin=1, norm=colors.LogNorm())
        ax2.set_xlabel("X")
        ax2.set_ylabel("Z")
        ax2.plot(om_x, om_z, 'x', color='r', markersize=1, alpha=0.2)
        
        c2 = plt.Circle((0, 0), 880, color='r', fill=False, linestyle='--')
        plt.gca().add_patch(c2)
        cbar2 = plt.colorbar(im1, ax=ax2)
        cbar2.set_label('Counts')
        ax2.set_aspect('equal')

        ax3 = fig.add_subplot(133, projection='3d')
        ax3.plot(first_mctree.col('x'),first_mctree.col('y'),first_mctree.col('z'), 
                 'o', markersize=0.2)#, alpha=0.3)
        
        fig.suptitle(f"{j}")

        plt.tight_layout()
        fig.savefig(f'./unweighted_plots/unweighted_injectpoint_{j}.png')
        del fig, ax1, ax2, ax3
                
        
    