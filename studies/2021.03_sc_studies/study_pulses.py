#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/icetray-start
#METAPROJECT combo/V01-01-00

import matplotlib.pyplot as plt
import numpy as np
import copy

from icecube import icetray, dataio, dataclasses
from icecube.icetray import I3Units
from I3Tray import I3Tray
import tools

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-i", type=str,
	dest="input_files",required=True,
	help="full path to the input file",
	)
parser.add_argument("-y", type=int,
	dest="year", default=2015,
	help="Which year of standard candle, e.g. 2015"
	)
parser.add_argument("-c", type=int,
	dest="candle", default=2,
	help="Which candle, e.g. 2"
	)
parser.add_argument("-f", type=int,
	dest="filter_setting", required=True,
	help="Which standard candle filter setting, e.g. 0, 1, 2, ..."
	)
# parser.add_argument("-o", type=str, 
# 	dest="output_dir",required=True,
# 	help="directory where the output should be written to"
# 	)
args = parser.parse_args()
filename = args.input_files
# output_location = args.output_dir

pulse_name='SplitInIcePulses'
# pulse_name='EHEBestPortiaPulseSRT'

frameId = 0
maxEvents = 10

scdata = '/misc/disk15/data/IceCube/RealData/86strings/standardcandle/2015/sc2/'
gcdfile = 'Level2_IC86.2014_data_Run00125920_0116_0_138_GCD.i3.gz'

geo = None

# om_strings = []
# om_depths = []

gcd_file_in = dataio.I3File(scdata + gcdfile)
while gcd_file_in.more():
	try:
		frame = gcd_file_in.pop_frame()
	except:
		continue
	
	if frame.Stop == icetray.I3Frame.Geometry:
		geo = frame['I3Geometry'].omgeo
		# for omkey, omgeo in geo:
		# 	if omgeo.position.z > 1000:
		# 		continue
		# 	om_strings.append(omkey.string)
		# 	om_depths.append(omgeo.position.z)

file_in = dataio.I3File(filename)
while file_in.more() and frameId < maxEvents:
	frameId+=1
	try:
		frame = file_in.pop_frame()
	except:
		continue

	if frame.Stop == icetray.I3Frame.DAQ:

		header = frame.Get('I3EventHeader')

		plot_waveforms = False
		if plot_waveforms:
			tools.print_waveforms(frame, "./plots",
				waveform_name='CalibratedWaveforms_ATWD',
				string=55, dom=44,
				title_mod=f'{args.year}_SC{args.candle}_F{args.filter_setting}'
				)
			tools.print_waveforms(frame, "./plots",
				waveform_name='CalibratedWaveforms_FADC',
				string=55, dom=44,
				title_mod=f'{args.year}_SC{args.candle}_F{args.filter_setting}'
				)

	if frame.Stop == icetray.I3Frame.Physics:
		header = frame.Get('I3EventHeader')
		if header.sub_event_stream!='InIceSplit':
			continue

		pulse_series = tools.get_pulse_series(frame, pulse_name)
		# print(pulse_series.keys())
		# exit(1)
		
		plot_panopticon = True
		if plot_panopticon:
			omkey_npe_dict = tools.get_homogqtot_omkey_npe_dict(pulse_series)

			plot_strings = []
			plot_depths = []
			plot_charges = []

			for omkey, omgeo in geo:
				if omgeo.position.z > 1000:
					continue
				plot_strings.append(omkey.string)
				plot_depths.append(geo[omkey].position.z)
				this_charge = -1
				if omkey in omkey_npe_dict:
					this_charge = omkey_npe_dict[omkey]
				plot_charges.append(this_charge)
				if omkey.string == 55 and omkey in omkey_npe_dict:
					print("OMKey {}, Depth {:.2f}".format(omkey, geo[omkey].position.z))

			plot_strings = np.asarray(plot_strings)
			plot_depths = np.asarray(plot_depths)
			plot_charges = np.asarray(plot_charges)

			hit_doms = plot_charges > 0


			# make a plot
			fig, axs = plt.subplots(1,1,figsize=(20,10))
			
			axs.scatter(plot_strings[~hit_doms], plot_depths[~hit_doms],
				marker='o', s=1, c='grey')

			sc2 = tools.get_sc_location(2)
			axs.plot(sc2['string'], sc2['z'], 'rx')

			this_ax = axs.scatter(plot_strings[hit_doms], plot_depths[hit_doms], 
				c=np.log10(plot_charges[hit_doms]))
			cbar = plt.colorbar(this_ax)
			cbar.set_label(r'$log_{10}$(NPE)')			

			axs.set_xlabel('String Number')
			axs.set_ylabel('Depth')
			axs.grid()
			# axs.set_ylim([-2800, -1200])
			axs.set_title(f'Ev {header.event_id}.{header.sub_event_id}, Year {args.year}, SC {args.candle}, Filter {args.filter_setting}')
			plt.tight_layout()
			fig.savefig('om_map_{}_{}_{}_{}.png'.format(header.run_id, header.event_id, frameId, pulse_name))
			print("------")

		plot_time_diff = False
		if plot_time_diff:
			launches = frame.Get('CleanInIceRawData')
			diffs = []
			for omkey, pulses in pulse_series.items():
				first_pulse_time = pulses[0].time
				first_launch_time = launches[omkey][0].time
				diffs.append(first_pulse_time - first_launch_time)
				# print("First pulse time {}, First launch time {}".format(first_pulse_time, first_launch_time))

			bins = np.arange(-200,0,5)
			fig, axs = plt.subplots(1,1,figsize=(5,5))
			axs.hist(diffs, bins=bins, alpha=0.5)
			axs.set_xlabel(r'$T_{firstpulse} - T_{launch}$')
			axs.set_ylabel('DOMs')
			plt.tight_layout()
			# axs.set_yscale('log')
			fig.savefig('time_diff_{}.png'.format(frameId))
			del fig, axs


