import h5py
import argparse
import numpy as np
import numpy
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import itertools

def find_contours_2D(x_values,y_values,xbins,weights=None,c1=16,c2=84):   
	"""
	Find upper and lower contours and median
	x_values = array, input for hist2d for x axis (typically truth)
	y_values = array, input for hist2d for y axis (typically reconstruction)
	xbins = values for the starting edge of the x bins (output from hist2d)
	c1 = percentage for lower contour bound (16% - 84% means a 68% band, so c1 = 16)
	c2 = percentage for upper contour bound (16% - 84% means a 68% band, so c2=84)
	Returns:
		x = values for xbins, repeated for plotting (i.e. [0,0,1,1,2,2,...]
		y_median = values for y value medians per bin, repeated for plotting (i.e. [40,40,20,20,50,50,...]
		y_lower = values for y value lower limits per bin, repeated for plotting (i.e. [30,30,10,10,20,20,...]
		y_upper = values for y value upper limits per bin, repeated for plotting (i.e. [50,50,40,40,60,60,...]
	"""
	if weights is not None:
		import wquantiles as wq
	y_values = numpy.array(y_values)
	indices = numpy.digitize(x_values,xbins)
	r1_save = []
	r2_save = []
	median_save = []
	for i in range(1,len(xbins)):
		mask = indices==i
		if len(y_values[mask])>0:
			if weights is None:
				r1, m, r2 = numpy.percentile(y_values[mask],[c1,50,c2])
			else:
				r1 = wq.quantile(y_values[mask],weights[mask],c1/100.)
				r2 = wq.quantile(y_values[mask],weights[mask],c2/100.)
				m = wq.median(y_values[mask],weights[mask])
		else:
			#print(i,'empty bin')
			r1 = numpy.nan
			m = numpy.nan
			r2 = numpy.nan
		median_save.append(m)
		r1_save.append(r1)
		r2_save.append(r2)
	median = numpy.array(median_save)
	lower = numpy.array(r1_save)
	upper = numpy.array(r2_save)

	x = list(itertools.chain(*zip(xbins[:-1],xbins[1:])))
	y_median = list(itertools.chain(*zip(median,median)))
	y_lower = list(itertools.chain(*zip(lower,lower)))
	y_upper = list(itertools.chain(*zip(upper,upper)))
	
	return x, y_median, y_lower, y_upper



parser = argparse.ArgumentParser()
parser.add_argument("-f", type=str, nargs='+',
	dest="input_files", required=True,
	help="paths to input files (absolute or relative)")
args = parser.parse_args()
files = args.input_files

reco_azimuth = np.asarray([])
true_azimuth = np.asarray([])

for file in files:

	file_in = h5py.File(file, "r")
	# print(file_in['LineFit'].dtype.names)

	try:
		reco_azimuth = np.concatenate((reco_azimuth, 
			file_in['LineFit']['azimuth']))
		true_azimuth = np.concatenate((true_azimuth,
			file_in['NuPrimary']['azimuth']))
		print("Finished {}".format(file))
	except:
		print('Skipping {}'.format(file))

	file_in.close()


reco_azimuth = np.rad2deg(reco_azimuth)
true_azimuth = np.rad2deg(true_azimuth)
bins = [np.linspace(0,360,180), np.linspace(0,360,180)]

# 2D histogram
fig = plt.figure(figsize=(6,5))
ax = fig.add_subplot(111)
counts, xedges, yedges, im = ax.hist2d(
	true_azimuth,
	reco_azimuth,
	bins=bins,
	cmin=1,
	norm=colors.LogNorm()
	)

# get the contours
x, y_med, y_lo, y_hi = find_contours_2D(
	x_values=true_azimuth,
	y_values=reco_azimuth,
	xbins=xedges
	)
# plot them
ax.plot(x, y_med, 'r-', label='Median')
ax.plot(x, y_lo, 'r-.', label='68% contour')
ax.plot(x, y_hi, 'r-.')

cbar = plt.colorbar(im, ax=ax)
cbar.set_label('Number of Events')#, fontsize=sizer)
ax.set_ylabel('Reco Azimuth')
ax.set_xlabel('True Azimuth')
ax.legend()
plt.tight_layout()
fig.savefig('test.png', dpi=300)
del fig, ax


