import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import style
import tables
import copy

style.use('/home/brian/IceCube/ehe/max_tools/EHE_analysis/eheanalysis/ehe.mplstyle')

# custom
from ehefluxes import fluxes
from eheanalysis import weighting, plotting
gzk_flux = fluxes.EHEFlux("ahlers_gzk")
from functools import partial
gzk_partial = partial(gzk_flux, which_species="nue_sum") 


def get_stuff(f):

    # get the charge
    charge = f.get_node('/Homogenized_QTot')
    charge = charge.col('value')

    # weight dict and prop matrix
    weight_dict = f.get_node('/JulietWeightDict')
    prop_matrix = [f.get_node(f'/PropagationMatrix{flav}') for flav in ['NuE', 'NuMu', 'NuTau']]

    evts_per_file = weight_dict.col('EventsPerFile')[0]

    return charge, weight_dict, prop_matrix, evts_per_file

species = ["nue", "numu", "nutau", "mu", "tau"]
# species = ["mu"]
file_paths = {
    "nue": "/home/mmeier/data/simulations/table_based_sim/juliet/nue/high_energy/l2_prime/nue_high_energy_l2_prime_1k_merged.hd5",
    "numu": "/home/mmeier/data/simulations/table_based_sim/juliet/numu/high_energy/l2_prime/numu_high_energy_l2_prime_1k_merged.hd5",
    "nutau": "/home/mmeier/data/simulations/table_based_sim/juliet/nutau/high_energy/l2_prime/nutau_high_energy_l2_prime_1k_merged.hd5",
    "mu": "/home/mmeier/data/simulations/table_based_sim/juliet/mu/high_energy/l2_prime/mu_high_energy_l2_prime_1k_merged.hd5",
    "tau": "/home/mmeier/data/simulations/table_based_sim/juliet/tau/high_energy/l2_prime/tau_high_energy_l2_prime_1k_merged.hd5"
}
# file_paths = {
#     "nue": "/home/mmeier/data/simulations/table_based_sim/juliet/nue/high_energy/l2_prime/1/Level2_prime_00000001.hd5",
#     "numu": "/home/mmeier/data/simulations/table_based_sim/juliet/numu/high_energy/l2_prime/1/Level2_prime_00000001.hd5",
#     "nutau": "/home/mmeier/data/simulations/table_based_sim/juliet/nutau/high_energy/l2_prime/1/Level2_prime_00000001.hd5",
#     "mu": "/home/mmeier/data/simulations/table_based_sim/juliet/mu/high_energy/l2_prime/1/Level2_prime_00000001.hd5",
#     "tau": "/home/mmeier/data/simulations/table_based_sim/juliet/tau/high_energy/l2_prime/1/Level2_prime_00000001.hd5",
# }

files = {
}

n_files = {
    "nue": 1000,
    "numu":  1000,
    "nutau":  1000,
    "mu":  1000,
    "tau":  1000,
}

livetime = 86400 * 365 # 1 year

charge_bins = np.logspace(2, 8, 51)
charge_bin_centers = (charge_bins[1:] + charge_bins[:-1]) * 0.5

# get weights by calculating them ourselves
def sample_flux(e_in_gev):
    return 1E-8 * (e_in_gev ** -2.)

enu_bins, enu_bin_centers = weighting.get_juliet_enu_binning()

summed_hist = None
hqtots = []
portias = []
e_dets = []

fig6, ax6 =  plt.subplots(1,1)
bins = np.linspace(4, 12, 60)

for s in species:

    print(f"Working on {s}")
    the_f = tables.open_file(file_paths[s])

    fname = file_paths[s]
    charge, weight_dict, prop_matrix, evts_per_file = get_stuff(the_f)
    n_gen = n_files[s] * evts_per_file

    for q in charge:
        hqtots.append(q)

    charge_portia = the_f.get_node("/EHEPortiaEventSummarySRT").col('bestNPEbtw')
    for q in charge_portia:
        portias.append(q)
    
    primary = the_f.get_node('/I3JulietPrimaryParticle')
    primary_e = primary.col('energy')
    ax6.hist(np.log10(primary_e), bins=bins, histtype='step', label=f'{s}: {len(primary_e)}')

    h = weighting.make_enu_2d_hist(
        weight_dict,
        n_gen,
        gzk_partial,
        var_values=charge,
        var_bins=charge_bins,
        prop_matrix=prop_matrix,
        livetime=livetime)
    
    if summed_hist is None:
        summed_hist = copy.deepcopy(h)
    else:
        summed_hist += copy.deepcopy(h)
    
    the_f.close() # close the file

ax6.set_yscale('log')
ax6.set_xlabel('Energy at Detector log10(GeV)')
ax6.set_ylabel('# Unweighted Events at L2')
ax6.legend()
plt.tight_layout()
fig6.savefig('edet.png')

hqtots = np.asarray(hqtots)
portias = np.asarray(portias)

# make projection along the charge
q_projection = h.sum(axis=1)
q_stepped_path = plotting.stepped_path(charge_bins, q_projection)

# calculate the efficiency
passing_rate = []
for q in charge_bin_centers:
    mask = charge_bin_centers > q
    passing_rate.append(np.sum(q_projection[mask]))
passing_rate = np.asarray(passing_rate)
passing_rate /= np.sum(q_projection)
the_q_cut = 27500
the_eff = np.interp(the_q_cut, charge_bin_centers, passing_rate)
print(f"The eff for {the_q_cut} is {the_eff}")

fig, ax = plt.subplots(1,1)
X, Y = np.meshgrid(enu_bin_centers, charge_bin_centers)
im = ax.pcolor(X, Y, h, norm=colors.LogNorm(vmin=1e-8))
ax.set_xscale('log')
ax.set_yscale('log')
cbar = plt.colorbar(im, ax=ax, label='Events / yr')
ax.set_xlabel(r'$E_{\nu}$ / GeV')
ax.set_ylabel('Q / PE')
fig.tight_layout()
fig.savefig('q_vs_e.png')

fig2, ax2 = plt.subplots(1,1)
ax2.plot(*q_stepped_path)
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_ylabel(r'Events / yr / bin')
ax2.set_xlabel('Q / PE')
ax2.set_ylim([1E-4, 2E-2])
fig2.tight_layout()
fig2.savefig('q_hist.png')

fig3, ax3 = plt.subplots(1,1)
ax3.plot(charge_bin_centers, passing_rate)
ax3.set_xscale('log')
ax3.set_ylabel(r'Efficiency')
ax3.set_xlabel('Q Cut / PE')
fig3.tight_layout()
fig3.savefig('q_cut_efficiency.png')

fig4, (ax4, ax5) = plt.subplots(1,2,figsize=(10,5))
# fig4, ax4 = plt.subplots(1,1)
bins = [np.linspace(2,8,51), np.linspace(2,8,51)]
vals, xedges, yedges, im = ax4.hist2d(
    np.log10(portias), np.log10(hqtots), bins=bins, 
    cmap=plt.cm.plasma,
    norm=colors.LogNorm()#, cmin=cmin
)
ax4.set_xlabel('log10(Portia Q)')
ax4.set_ylabel('log10(HQtot)')
ax4.plot([2,8],[2,8], '--')
charge_cbar = plt.colorbar(im, ax=ax4, label="Counts")
im.set_clim(1,1E5)

cmap = plt.cm.plasma
norm = colors.LogNorm()
bins = [np.linspace(2,8,51), np.linspace(-1,1,51)]
vals, xedges, yedges, im = plotting.make_2D_hist(
    ax5, np.log10(portias), hqtots/portias-1., bins,
    cmap=cmap, norm=norm, cmin=1.,
    xlabel='log10(Portia Q)', ylabel='HQtot/Portia - 1'
)
x, y_med, y_lo, y_hi = plotting.find_contours_2D(
    x_values=np.log10(portias),
    y_values=hqtots/portias-1.,
    xbins=xedges, c1=0, c2=68
)
ax5.plot(x, y_med, 'r-', label='Median')
ax5.plot(x, y_hi, 'r-.', label='68% contour')
ax5.plot(x, y_lo, 'r-.')
ax5.plot([2,8],[0,0], '--')
ratio_cbar = plt.colorbar(im, ax=ax5, label="Counts")
im.set_clim(1,1E5)

fig4.tight_layout()
fig4.savefig('hqtot_vs_portia.png')