import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import argparse

import utils_weights
livetime = 86400 * 365


def get_G_frame_stuff(gcd_file):
	from icecube import dataio, dataclasses, icetray
	geo = None
	file = dataio.I3File(str(gcd_file))
	while(file.more()):
		frame = file.pop_frame()
		if frame.Stop == icetray.I3Frame.Geometry:
			geo = frame['I3Geometry'].omgeo
		else:
			continue
	if geo is None:
		raise LookupError('Cannot find a G frame in the provided file {}'.format(gcd_file))

	om_x = []
	om_z = []
	strings_x = []
	strings_y = []
	finished = []
	for omkey, omgeo in geo:
		if omgeo.position.z > 1000:
			continue
		om_x.append(omgeo.position.x)
		om_z.append(omgeo.position.z)
		if omkey.om == 1:
			if omkey.string not in finished:
				finished.append(omkey.string)
				strings_x.append(omgeo.position.x)
				strings_y.append(omgeo.position.y)
	om_x = np.asarray(om_x)
	om_z = np.asarray(om_z)
	strings_x = np.asarray(strings_x)
	strings_y = np.asarray(strings_y)
	z_surface = dataclasses.I3Constants.SurfaceElev - dataclasses.I3Constants.OriginElev
	# om_z -= z_surface
	return om_x, om_z, strings_x, strings_y, z_surface

parser = argparse.ArgumentParser()
parser.add_argument("-f", type=str,
	dest="file", required=True,
	help="paths to numu file")
parser.add_argument("-g", type=str,
    dest="gcd_file",
    default='/cvmfs/icecube.opensciencegrid.org/data/GCD/GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE_withScaledNoise.i3.gz',
    help='gcd file',
    )
args = parser.parse_args()

om_x, om_z, strings_x, strings_y, surface = get_G_frame_stuff(args.gcd_file)

atmo_flux_model = utils_weights.get_flux_model('H3a_SIBYLL23C', 'nugen')
cr_flux_model = utils_weights.get_flux_model('GaisserH3a', 'corsika')

file = pd.HDFStore(args.file)

numu_weighter = utils_weights.get_weighter(file, 'nugen', 1000)
numu_atmo_weights = numu_weighter.get_weights(atmo_flux_model)
numu_xx = numu_weighter.get_column('VertexPosition', 'x')
numu_yy = numu_weighter.get_column('VertexPosition', 'y')
numu_zz = numu_weighter.get_column('VertexPosition', 'z')
file.close()

numu_atmo_weights *= livetime

cmap=plt.cm.viridis
sizer=15

do_event_dist=True
if do_event_dist:

	scale=1000

	fig = plt.figure(figsize=(18,7))
	bins = [np.linspace(-20000/scale,20000/scale,4000), np.linspace(-20000/scale,20000/scale,4000)]
	bins = 100
	ax = fig.add_subplot(121)
	counts, xedges, yedges, im = ax.hist2d(
			numu_xx/scale, 
			numu_yy/scale, 
			bins=bins,
			cmap=cmap,
			norm=colors.LogNorm(),
			cmin=1
			)
	cbar = plt.colorbar(im, ax=ax)
	ax.plot(strings_x/scale, strings_y/scale, 'x', markersize=1, color='r')
	ax.set_xlim([-10000/scale, 10000/scale])
	ax.set_ylim([-10000/scale, 10000/scale])
	cbar.set_label('Events', fontsize=sizer)
	ax.set_ylabel(r'Y [km]', fontsize=sizer)
	ax.set_xlabel(r'X [km]', fontsize=sizer)
	ax.tick_params(labelsize=sizer)
	ax.set_aspect('equal')
	ax.set_title('Top-Down View')

	bins = [np.linspace(-20000/scale,20000/scale,400), np.linspace(-20000/scale,20000/scale,400)]
	ax2 = fig.add_subplot(122)
	counts, xedges, yedges, im = ax2.hist2d(
			numu_xx/scale, 
			numu_zz/scale, 
			bins=bins,
			cmap=cmap,
			norm=colors.LogNorm(),
			cmin=1
			)
	cbar2 = plt.colorbar(im, ax=ax2)
	ax2.plot(om_x/scale, om_z/scale, 'x', markersize=1, color='r')
	ax2.axhline(y=surface/scale, xmin=-20000/scale, xmax=20000/scale, color='r', linestyle='--')
	ax2.set_xlim([-10000/scale, 10000/scale])
	ax2.set_ylim([-5000/scale, 2000/scale])
	# cbar2.set_label('Events', fontsize=sizer)
	ax2.set_xlabel(r'X [km]', fontsize=sizer)
	ax2.set_ylabel(r'Z [km]', fontsize=sizer)
	ax2.tick_params(labelsize=sizer)
	ax2.set_title('Side view')
	ax2.set_aspect('equal')


	plt.tight_layout()
	fig.savefig('event_dist.png', edgecolor='none', bbox_inches='tight', dpi=300)
