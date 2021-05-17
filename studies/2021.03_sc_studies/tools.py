from icecube import icetray, dataio, dataclasses, recclasses
from icecube.icetray import OMKey, I3Units
import matplotlib.pyplot as plt
import numpy as np
import os

# return the MJD of the time windows
# these times were chosen by carefully selecting where the cut should go 
# to exclude the times when the filter wheel rotates to 100%
def get_start_stop(year, candle, filter):
	offset= 0
	start = 0
	stop = 0
	if year==2015:
		if candle==2:
			offset = 5.7038e4
			if filter==1:
				start = 0.0645
				stop = 0.0878
				charge_min = 0.5E5
				charge_max = 0.8E5
			if filter==2:
				start = 0.08855
				stop = 0.113
				charge_min = 0.8E5
				charge_max = 1.2E5
			if filter==3:
				start = 0.114
				stop = 0.139
				charge_min = 1.5E5
				charge_max = 2.2E5
			if filter==4:
				start = 0.14
				stop = 0.165
				charge_min = 2.8E5
				charge_max = 4E5
			if filter==5:
				start = 0.1678
				stop = 0.191
				charge_min = 3.6E5
				charge_max = 4.7E5
			if filter==6:
				start = 0.192
				stop = 0.22
				charge_min = 6.1E5
				charge_max = 7.5E5
	return offset+start, offset+stop, charge_min, charge_max

'''
	Decide where to cut an event based on it's time and charge
	start and stop are I3Times
	qmin and qmax are homogenized Q tot variables in PE
'''

def cut_by_config(frame, start, stop, qmin, qmax):

	keeper = False
	if frame.Has('I3EventHeader'):
		evt_time = frame.Get('I3EventHeader').start_time
		if evt_time > start and evt_time < stop:
			keeper = True
			# don't impose charge cut any more
			# if frame.Has('HomogenizedQTot'):
				# hqtot = frame.Get('HomogenizedQTot').value
				# if hqtot > qmin and hqtot < qmax:
					# print("keep, qtot {}".format(hqtot))
					# keeper = True

	return keeper

def get_pulse_series(frame,pulse_mask_name):
	hit_map = None
	if type(frame[pulse_mask_name]) == dataclasses.I3RecoPulseSeriesMap:
		hit_map = frame.Get(pulse_mask_name)
	elif type(frame[pulse_mask_name]) == dataclasses.I3RecoPulseSeriesMapMask:
		hit_map = frame[pulse_mask_name].apply(frame)
	else:
		# just get something, the user better deal with it...
		hit_map = frame.Get(pulse_mask_name)
	return hit_map


def get_homogqtot_omkey_npe_dict(pulse_series):
	omkey_npe_dict = {}
	for omkey, pulses in pulse_series:
		charge_this_dom = 0
		if type(pulses) == recclasses.I3PortiaPulse:
			pulses = [pulses] # put this into a list, just in case it's a PortiaPulse
		for p in pulses:
			if type(p) == dataclasses.I3RecoPulse:
				charge_this_dom += p.charge
			elif type(p) == recclasses.I3PortiaPulse:
				charge_this_dom += p.GetEstimatedNPE()
		omkey_npe_dict[omkey] = charge_this_dom
	return omkey_npe_dict

# see https://github.com/icecube/icetray/blob/main/clsim/python/StandardCandleFlasherPulseSeriesGenerator.py#L75
def get_sc_location(candle):
	if candle==1:
		return {"x": 544.07, "y": 55.89,"z": 136.86, "string": 40}
	elif candle==2:
		return {"x": 11.87, "y": 179.19,"z": -205.64, "string": 55}

def print_waveforms(frame, outputdir, waveform_name, string, dom,
	title_mod=None):

	header = frame['I3EventHeader']
	key = OMKey(string, dom)
	try:
		waveform_series = frame.Get(waveform_name)[key]
	except:
		print('no waveforms found for this dom. skip')
		return

	# make plots
	fig, ax = plt.subplots()
	ax.set_xlabel(r'Time / ns')

	min_time = 0
	for waveform in waveform_series:
		wf_vect = np.array(waveform.waveform) / I3Units.mV
		start_time = waveform.time
		bin_width = waveform.bin_width
		if min_time==0:
			min_time = start_time

		# Skip possible second launches
		if len(wf_vect) == 128:
			if start_time > min_time + 5000:
				continue
		elif len(wf_vect) == 256:
			if start_time > min_time + 20000:
				continue

		time = np.linspace(start_time, start_time + bin_width * len(wf_vect), len(wf_vect))

		ax.plot(time, wf_vect, label=f'{waveform.source}, Ch {waveform.channel}')
		if len(wf_vect) == 128:
			ax.set_ylabel(r'Voltage / mV')
		elif len(wf_vect) == 256:
			ax.set_ylabel('Voltage / mV')

	ax.set_title(f'({string}, {dom}), Ev {header.event_id}.{header.sub_event_id}, {waveform_name} ({title_mod})')

	ax.legend(loc='best')
	if len(wf_vect) == 128:
		ax.set_ylim([-50,7000])
	filename = os.path.join(
		outputdir,
		f'run{header.run_id}_evt{header.event_id}_sevt_{header.sub_event_id}_str{string}_dom{dom}_{title_mod}.png')
	if len(wf_vect) == 256:
		filename = filename.replace('.png', '_fadc.png')
	print(filename)
	fig.savefig(filename)



