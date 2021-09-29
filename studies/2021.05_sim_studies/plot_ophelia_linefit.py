import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import argparse, os

import utils_weights
import utils_plot
livetime = 86400 * 365

parser = argparse.ArgumentParser()
parser.add_argument("-file", type=str,
    dest="file", required=True,
    help="paths to numu file")
args = parser.parse_args()
file_name = args.file
file_name = os.path.basename(file_name)
dataset = file_name.split('_')[1].split('.')[0]

in_file = pd.HDFStore(args.file)

atmo_flux_model = utils_weights.get_flux_model('H3a_SIBYLL23C', 'nugen')

weighter = utils_weights.get_weighter(in_file, 'nugen', 1000)

ophelia_zenith = weighter.get_column('EHEOpheliaParticleSRT_ImpLF', 'zenith')
ophelia_azimuth = weighter.get_column('EHEOpheliaParticleSRT_ImpLF', 'azimuth')
linefit_zenith = weighter.get_column('LineFit', 'zenith')
linefit_azimuth = weighter.get_column('LineFit', 'azimuth')
true_zenith = weighter.get_column('PolyplopiaPrimary', 'zenith')
true_azimuth = weighter.get_column('PolyplopiaPrimary', 'azimuth')
energy = weighter.get_column('PolyplopiaPrimary', 'energy')
hqtot = weighter.get_column('Homogenized_QTot', 'value')
weights = weighter.get_weights(atmo_flux_model)

in_file.close()
weights *= livetime

# cmap=plt.cm.plasma
cmap=plt.cm.viridis
sizer=15
kwargs = {'cmap': cmap, 
            'norm' : colors.LogNorm(), 
            'cmin': 1,
            # 'weights' : weights
            }
kwargs_zenith = kwargs.copy()
bins =  [np.linspace(0,1,100), np.linspace(0,1,100)]
kwargs_zenith['bins'] = bins
if int(dataset) == 21218:
    bins_e = np.linspace(3,8,40)
else:
    bins_e = np.linspace(4,9,40)
mask = np.log10(hqtot) > 3

def cartesian_components(zenith, azimuth):
    theta = zenith
    phi = azimuth
    return -np.sin(theta)*np.cos(phi), -np.sin(theta)*np.sin(phi), -np.cos(theta)

# broadcasting was numpy's greatest idea
ophelia_x, ophelia_y, ophelia_z = cartesian_components(ophelia_zenith, ophelia_azimuth)
linefit_x, linefit_y, linefit_z = cartesian_components(linefit_zenith, linefit_azimuth)
true_x, true_y, true_z = cartesian_components(true_zenith, true_azimuth)
opening_angle_ophelia = np.degrees(np.arccos(ophelia_x*true_x + ophelia_y*true_y + ophelia_z*true_z))
opening_angle_linefit = np.degrees(np.arccos(linefit_x*true_x + linefit_y*true_y + linefit_z*true_z))

opening_angle_ophelia = opening_angle_ophelia[mask]
opening_angle_linefit = opening_angle_linefit[mask]
ophelia_zenith = ophelia_zenith[mask]
linefit_zenith = linefit_zenith[mask]
true_zenith = true_zenith[mask]

do_energy_bins=True
if do_energy_bins:

    # energy distribution

    # bins_e = np.geomspace(1e5, 1e10, 50)
    fig = plt.figure(figsize=(4,4))
    ax = fig.add_subplot(111)
    ax.hist(np.log10(energy), bins=bins_e, alpha=0.5, label='All Events')
    ax.hist(np.log10(energy[mask]), bins=bins_e, alpha=0.5, label='HQTot > 1k')
    ax.set_xlabel(r'True Neutrino Energy log$_{10}$(GeV)')
    ax.set_ylabel(r'Events')
    ax.set_yscale('log')
    ax.legend()
    ax.set_title('{}'.format(dataset))
    plt.tight_layout()
    fig.savefig('plots/{}_energy_dist.png'.format(dataset), 
            edgecolor='none', bbox_inches='tight', dpi=300)
    del fig, ax


do_opening_angle=False
if do_opening_angle:

    # fig = plt.figure(figsize=(7,7))
    # ax = fig.add_subplot(111)
    # # bins = np.linspace(0,180,180)
    # bins_oa=100
    # ax.hist(opening_angle_ophelia, bins=bins_oa, alpha=0.5, label='Ophelia')
    # ax.hist(opening_angle_linefit, bins=bins_oa, alpha=0.5, label='LineFit')
    # ax.set_xlabel(r'Opening Angle [deg]')
    # ax.set_ylabel(r'Events')
    # ax.set_yscale('log')
    # ax.set_xscale('log')
    # ax.legend()
    # plt.tight_layout()
    # fig.savefig('plots/{}_opening_angle.png'.format(dataset), 
    #         edgecolor='none', bbox_inches='tight', dpi=300)
    # del fig, ax


    #------------------------------------
    #------------------------------------
    #------------------------------------

    # opening angle vs energy, 2D

    fig = plt.figure(figsize=(15,7))    
    ax = fig.add_subplot(121)
    counts, xedges, yedges, im = ax.hist2d(
            np.log10(energy[mask]),
            opening_angle_ophelia,
            bins=[bins_e, np.linspace(0,180,180)],
            **kwargs
            )
    cbar = plt.colorbar(im, ax=ax)
    x, y_med, y_lo, y_hi = utils_plot.find_contours_2D(
        x_values=np.log10(energy[mask]),
        y_values=opening_angle_ophelia,
        xbins=xedges,
        c1=0, c2=68,
        )
    ax.plot(x, y_med, 'r-', label='Median')
    ax.plot(x, y_hi, 'r-.', label='68% contour')
    ax.set_title('Ophelia ({})'.format(dataset))
    cbar.set_label('Events', fontsize=sizer)
    ax.set_ylabel(r'Ophelia Opening Angle', fontsize=sizer)
    ax.set_xlabel(r'Energy log$_{10}$(GeV)', fontsize=sizer)
    ax.tick_params(labelsize=sizer)
    ax.legend()

    ax2 = fig.add_subplot(122)
    counts, xedges, yedges, im = ax2.hist2d(
            np.log10(energy[mask]),
            opening_angle_linefit,
            bins=[bins_e, np.linspace(0,180,180)],
            **kwargs
            )
    cbar2 = plt.colorbar(im, ax=ax2)
    x, y_med, y_lo, y_hi = utils_plot.find_contours_2D(
        x_values=np.log10(energy[mask]),
        y_values=opening_angle_linefit,
        xbins=xedges,
        c1=0, c2=68,
        )
    ax2.plot(x, y_med, 'r-', label='Median')
    ax2.plot(x, y_hi, 'r-.', label='68% contour')
    ax2.set_title('LineFit ({})'.format(dataset))
    cbar2.set_label('Events', fontsize=sizer)
    ax2.set_ylabel(r'LineFit Opening Angle', fontsize=sizer)
    ax2.set_xlabel(r'Energy log$_{10}$(GeV)', fontsize=sizer)
    ax2.tick_params(labelsize=sizer)
    ax2.legend()

    plt.tight_layout()
    fig.savefig('plots/{}_opening_angle_vs_e_2d.png'.format(dataset), 
            edgecolor='none', bbox_inches='tight', dpi=300)
    del fig, ax, ax2

    #------------------------------------
    #------------------------------------
    #------------------------------------

    # opening angle vs energy, 1D
    
    fig = plt.figure(figsize=(7,4))    
    ax = fig.add_subplot(111)
    x, y_med, y_lo, y_hi = utils_plot.find_contours_2D(
        x_values=np.log10(energy[mask]),
        y_values=opening_angle_ophelia,
        xbins=xedges,
        c1=0, c2=68,
        )
    ax.plot(x, y_med, 'C0o-', label='Ophelia Median')
    ax.plot(x, y_hi, 'C0o--', label='Ophelia 68%')
    x, y_med, y_lo, y_hi = utils_plot.find_contours_2D(
        x_values=np.log10(energy[mask]),
        y_values=opening_angle_linefit,
        xbins=xedges,
        c1=0, c2=68
        )
    sizer=12
    ax.plot(x, y_med, 'C1s-', label='LineFit Median')
    ax.plot(x, y_hi, 'C1s--', label='LineFit 68%')
    ax.set_ylabel(r'$\Delta \Psi$ (Opening Angle)', fontsize=sizer)
    ax.set_xlabel(r'Energy log$_{10}$(GeV)', fontsize=sizer)
    ax.tick_params(labelsize=sizer)
    ax.legend()
    ax.set_ylim([0,40])
    if dataset == '21218':
        ax.set_ylim([60,100])
    ax.set_title('{}'.format(dataset))
    plt.tight_layout()
    fig.savefig('plots/{}_opening_angle_vs_e.png'.format(dataset), 
            edgecolor='none', bbox_inches='tight', dpi=300)
    del fig, ax
    sizer=15

do_zenith = True
if do_zenith:

    #------------------------------------
    #------------------------------------
    #------------------------------------

    # true vs reco cos(zen), 2D


    fig2 = plt.figure(figsize=(7,7))
    f2_ax = fig2.add_subplot(111)

    fig = plt.figure(figsize=(22,7))    

    # linefit vs ophelia
    ax = fig.add_subplot(131)
    counts, xedges, yedges, im = ax.hist2d(
            np.cos(ophelia_zenith),
            np.cos(linefit_zenith),
            **kwargs_zenith
            )
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Events', fontsize=sizer)
    ax.set_ylabel(r'LineFit $\cos(\theta)$', fontsize=sizer)
    ax.set_xlabel(r'Ophelia $\cos(\theta)$', fontsize=sizer)
    ax.tick_params(labelsize=sizer)
    ax.set_title('LineFit vs Ophelia ({})'.format(dataset))
    ax.plot([-1,1], [-1, 1], 'w:')
    x, y_med, y_lo, y_hi = utils_plot.find_contours_2D(
        x_values=np.cos(ophelia_zenith),
        y_values=np.cos(linefit_zenith),
        xbins=xedges
        )
    ax.plot(x, y_med, 'r-', label='Median')
    ax.plot(x, y_lo, 'r-.', label='68% contour')
    ax.plot(x, y_hi, 'r-.')
    ax.legend()
    ax.set_aspect('equal')

    # ophelia vs truth
    ax2 = fig.add_subplot(132)
    counts, xedges, yedges, im = ax2.hist2d(
            np.cos(true_zenith),
            np.cos(ophelia_zenith),
            **kwargs_zenith
            )
    cbar2 = plt.colorbar(im, ax=ax2)
    cbar2.set_label('Events', fontsize=sizer)
    ax2.set_xlabel(r'True $\cos(\theta)$', fontsize=sizer)
    ax2.set_ylabel(r'Ophelia $\cos(\theta)$', fontsize=sizer)
    ax2.tick_params(labelsize=sizer)
    ax2.set_title('Ophelia vs Truth ({})'.format(dataset))
    ax2.plot([-1,1], [-1, 1], 'w:')
    x, y_med, y_lo, y_hi = utils_plot.find_contours_2D(
        x_values=np.cos(true_zenith),
        y_values=np.cos(ophelia_zenith),
        xbins=xedges
        )
    ax2.plot(x, y_med, 'r-', label='Median')
    ax2.plot(x, y_lo, 'r-.', label='68% contour')
    ax2.plot(x, y_hi, 'r-.')
    ax2.set_aspect('equal')
    f2_ax.plot(x, y_med - x, 'o', label='Ophelia')

    # linefit vs truth
    ax3 = fig.add_subplot(133)
    counts, xedges, yedges, im = ax3.hist2d(
            np.cos(true_zenith),
            np.cos(linefit_zenith),
            **kwargs_zenith
            )
    cbar3 = plt.colorbar(im, ax=ax3)
    cbar3.set_label('Events', fontsize=sizer)
    ax3.set_xlabel(r'True $\cos(\theta)$', fontsize=sizer)
    ax3.set_ylabel(r'LineFit $\cos(\theta)$', fontsize=sizer)
    ax3.tick_params(labelsize=sizer)
    ax3.set_title('LineFit vs Truth ({})'.format(dataset))
    ax3.plot([-1,1], [-1, 1], 'w:')
    x, y_med, y_lo, y_hi = utils_plot.find_contours_2D(
        x_values=np.cos(true_zenith),
        y_values=np.cos(linefit_zenith),
        xbins=xedges
        )
    ax3.plot(x, y_med, 'r-', label='Median')
    ax3.plot(x, y_lo, 'r-.', label='68% contour')
    ax3.plot(x, y_hi, 'r-.')
    ax3.set_aspect('equal')
    f2_ax.plot(x, y_med - x, 's', label='LineFit')

    plt.tight_layout()
    fig.savefig('plots/{}_ophelia_vs_linefit.png'.format(dataset), 
        edgecolor='none', bbox_inches='tight', dpi=300)


    f2_ax.set_xlabel(r'True $\cos(\theta)$')
    f2_ax.set_ylabel(r'Median Error in $\cos(\theta)$')
    f2_ax.set_ylim([-0.2, 0.2])
    if dataset == '21218':
        f2_ax.set_ylim([-1, 1])
    f2_ax.set_title('{}'.format(dataset))
    f2_ax.legend()
    fig2.savefig('plots/{}_ophelia_vs_linefit_med_error.png'.format(dataset),
        edgecolor='none', bbox_inches='tight', dpi=300
        )
    del fig, ax, fig2, f2_ax

    #------------------------------------
    #------------------------------------
    #------------------------------------

    # cos(zen) resolution vs energy

    fig = plt.figure(figsize=(15,7))    
    ax = fig.add_subplot(121)
    counts, xedges, yedges, im = ax.hist2d(
            np.log10(energy[mask]),
            np.cos(true_zenith) - np.cos(ophelia_zenith),
            bins=[bins_e, np.linspace(-2,2,100)],
            **kwargs
            )
    cbar = plt.colorbar(im, ax=ax)
    x, y_med, y_lo, y_hi = utils_plot.find_contours_2D(
        x_values=np.log10(energy[mask]),
        y_values=np.cos(true_zenith) - np.cos(ophelia_zenith),
        xbins=xedges,
        # c1=0, c2=68,
        )
    ax.plot(x, y_med, 'r-', label='Median')
    ax.plot(x, y_hi, 'r-.', label='68% contour')
    ax.plot(x, y_lo, 'r-.')
    ax.set_title('Ophelia ({})'.format(dataset))
    cbar.set_label('Events', fontsize=sizer)
    ax.set_ylabel(r'$\Delta (\cos \theta) (True - Reco)$', fontsize=sizer)
    ax.set_xlabel(r'Energy log$_{10}$(GeV)', fontsize=sizer)
    ax.tick_params(labelsize=sizer)
    ax.legend()

    ax2 = fig.add_subplot(122)
    counts, xedges, yedges, im = ax2.hist2d(
            np.log10(energy[mask]),
            np.cos(true_zenith) - np.cos(linefit_zenith),
            bins=[bins_e, np.linspace(-2,2,100)],
            **kwargs
            )
    cbar2 = plt.colorbar(im, ax=ax2)
    x, y_med, y_lo, y_hi = utils_plot.find_contours_2D(
        x_values=np.log10(energy[mask]),
        y_values=np.cos(true_zenith) - np.cos(linefit_zenith),
        xbins=xedges,
        # c1=0, c2=68,
        )
    ax2.plot(x, y_med, 'r-', label='Median')
    ax2.plot(x, y_hi, 'r-.', label='68% contour')
    ax2.plot(x, y_lo, 'r-.')
    ax2.set_title('LineFit ({})'.format(dataset))
    cbar2.set_label('Events', fontsize=sizer)
    ax2.set_ylabel(r'$\Delta (\cos \theta) (True - Reco)$', fontsize=sizer)
    ax2.set_xlabel(r'Energy log$_{10}$(GeV)', fontsize=sizer)
    ax2.tick_params(labelsize=sizer)
    ax2.legend()

    plt.tight_layout()
    fig.savefig('plots/{}_deltazenith_vs_e_2d.png'.format(dataset), 
            edgecolor='none', bbox_inches='tight', dpi=300)
    del fig, ax, ax2

    #------------------------------------
    #------------------------------------
    #------------------------------------

    # cos(zen) resolution vs energy, 1D


    fig = plt.figure(figsize=(7,4))    
    ax = fig.add_subplot(111)
    x, y_med, y_lo, y_hi = utils_plot.find_contours_2D(
        x_values=np.log10(energy[mask]),
        y_values=np.cos(true_zenith) - np.cos(ophelia_zenith),
        xbins=xedges,
        # c1=0, c2=68,
        )
    # ax.plot(x, y_med, 'C0o-', label='Ophelia Median')
    # ax.plot(x, y_hi, 'C0o--', label='Ophelia 68%')
    # ax.errorbar(x, y_med, yerr=[y_hi-y_med, y_med-y_lo], capsize=0.0, fmt='o', label='Ophelia')
    ax.plot(x,y_med, color='C0', alpha=1, label='Ophelia Median')
    ax.fill_between(x,y_med, y_hi,color='C0', alpha=0.2, label='Ophelia 68%')
    ax.fill_between(x,y_med, y_lo, color='C0', alpha=0.2)

    x, y_med, y_lo, y_hi = utils_plot.find_contours_2D(
        x_values=np.log10(energy[mask]),
        y_values=np.cos(true_zenith) - np.cos(linefit_zenith),
        xbins=xedges,
        # c1=0, c2=68
        )
    sizer=12
    # ax.plot(x, y_med, 'C1s-', label='LineFit Median')
    # ax.plot(x, y_hi, 'C1s--', label='LineFit 68%')
    # ax.errorbar(x, y_med, yerr=[y_hi-y_med, y_med-y_lo], capsize=0.0, fmt='o', label='LineFit')
    ax.plot(x,y_med, color='C1', alpha=1, label='LineFit Median')
    ax.fill_between(x,y_med, y_hi,color='C1', alpha=0.2, label='LineFit 68%')
    ax.fill_between(x,y_med, y_lo, color='C1', alpha=0.2)

    ax.set_ylabel(r'$\Delta (\cos \theta) (True - Reco)$', fontsize=sizer)
    ax.set_xlabel(r'Energy log$_{10}$(GeV)', fontsize=sizer)
    ax.tick_params(labelsize=sizer)
    ax.legend(loc='upper left')
    ax.set_ylim([-0.25,1])
    if dataset == '21218':
        ax.set_ylim([-1, 1])
    ax.set_title('{}'.format(dataset))
    
    ax.plot([3,9], [0, 0], 'k:')
    plt.tight_layout()
    fig.savefig('plots/{}_deltazenith_angle_vs_e.png'.format(dataset), 
            edgecolor='none', bbox_inches='tight', dpi=300)
    del fig, ax
    sizer=15

