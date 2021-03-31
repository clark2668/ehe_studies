from icecube import icetray, dataio, dataclasses

# return the MJD of the time windows
def get_start_stop(year, candle, filter):
	offset= 0
	start = 0
	stop = 0
	if year==2015:
		if candle==2:
			offset = 5.7038e4
			if filter==1:
				start = 0.060
				stop = 0.090
				charge_min = 0.5E5
				charge_max = 0.8E5
			if filter==2:
				start = 0.085
				stop = 0.115
				charge_min = 0.8E5
				charge_max = 1.2E5
			if filter==3:
				start = 0.11
				stop = 0.14
				charge_min = 1.5E5
				charge_max = 2.2E5
			if filter==4:
				start = 0.1375
				stop = 0.165
				charge_min = 2.8E5
				charge_max = 4E5
			if filter==5:
				start = 0.166
				stop = 0.195
				charge_min = 3.6E5
				charge_max = 4.7E5
			if filter==6:
				start = 0.190
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
			if frame.Has('HomogenizedQTot'):
				hqtot = frame.Get('HomogenizedQTot').value
				if hqtot > qmin and hqtot < qmax:
					print("keep, qtot {}".format(hqtot))
					keeper = True

	return keeper

def get_pulse_series(frame,pulse_mask_name):
	if type(frame[pulse_mask_name]) == dataclasses.I3RecoPulseSeriesMap:
		hit_map = frame.Get(pulse_mask_name)
	elif type(frame[pulse_mask_name]) == dataclasses.I3RecoPulseSeriesMapMask:
		hit_map = frame[pulse_mask_name].apply(frame)
	return hit_map

def get_homogqtot_omkey_npe_dict(pulse_series):
	omkey_npe_dict = {}
	for omkey, pulses in pulse_series:
		charge_this_dom = 0
		for p in pulses:
			charge_this_dom += p.charge
		omkey_npe_dict[omkey] = charge_this_dom
	return omkey_npe_dict