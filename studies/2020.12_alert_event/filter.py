from icecube import icetray, dataclasses, portia
from icecube.icetray import I3Units
from icecube.phys_services.which_split import which_split
from I3Tray import Inf

import numpy as np
import h5py as h5py

# useful for helping to define the time of arrival of the biggest pulse
def find_time_of_largest_pulse(portia_pulse_map):
	biggest_pulse = 0.
	biggeset_pulse_time = 0.
	for omkey, pulses in portia_pulse_map:
		this_npe = pulses.GetEstimatedNPE()
		this_t10 = pulses.GetRecoPulse().time
		if this_npe > biggest_pulse:
			biggest_pulse = this_npe
			biggeset_pulse_time = this_t10
	return biggeset_pulse_time

# get pulse value, and also allow for window cleaning
def get_portia_pulse_values(portia_pulse, largest_time, window_start, window_end, enforce_window):
	npe = portia_pulse.GetEstimatedNPE()
	t10 = portia_pulse.GetRecoPulse().time
	npe_out = 0
	if enforce_window and ((t10 - largest_time >= window_start) and (t10 - largest_time <= window_end)):
		npe_out = npe
	else:
		npe_out = 0
	return npe_out

def get_portia_omkey_baseline_dict(splitted_dom_map, pulse_map):
	'''
	Get the baseline values for each omkey
	splitted_dom_map is the launched omkeys
	the pulse_map is a map of Portia pulses--can be fadc, atwd, or "best"
	'''

	omkey_baseline_dict = {}
	for omkey, launches in splitted_dom_map:

		if omkey in pulse_map:
			pulse = pulse_map[omkey]
			baseline = pulse.GetBaseLine()
			# print("Omkey {}, baseline {}".format(omkey, baseline))
			omkey_baseline_dict[omkey] = baseline

	return omkey_baseline_dict

def get_portia_omkey_npe_dict(splitted_dom_map, fadc_pulse_map, atwd_pulse_map, 
	high_qe_doms = [], doBTW=True, excludeFADC=False, excludeATWD=False):
	'''
	This is a roughly replication of the MakePortiaEvent function in I3Portia.xx
	We check all of the launched omkeys in the splitted_dom_map
	And get the FADC and ATWD pulse values for that omkey
	And first figure out if they are inside a basetime window (BTW).
	If they are, we figure out which is larger--the FADC or the ATWD--and use that
	as the NPE estimate for that event.
	We will return a dict of the omkey to the NPE estimate (the omkey_npe_dict)
	'''
	if excludeFADC and excludeATWD:
		icetray.logging.log_fatal("excludeFADC={} and excludeATWD={}, which doesn't make sense. Abort!".format(excludeFADC, excludeATWD))

	# now set up window for cleaning
	largest_time_fadc = find_time_of_largest_pulse(fadc_pulse_map)
	largest_time_best = find_time_of_largest_pulse(atwd_pulse_map)
	largest_time = max(largest_time_fadc, largest_time_best)
	start_time_btw = -4400.0*I3Units.ns
	end_time_btw = 6400.0*I3Units.ns
	
	best_npe = 0.
	omkey_npe_dict = {}

	for omkey, launches in splitted_dom_map:

		# skip the high qe doms
		if omkey in high_qe_doms:
			continue
		
		# fadc estimate for this event
		this_fadc = 0
		if omkey in fadc_pulse_map:
			fadc_pulse = fadc_pulse_map[omkey]
			this_fadc = get_portia_pulse_values(fadc_pulse, largest_time, 
				start_time_btw, end_time_btw, doBTW)

		# "best" estimate for this event
		this_atwd = 0
		if omkey in atwd_pulse_map:
			atwd_pulse = atwd_pulse_map[omkey]
			this_atwd = get_portia_pulse_values(atwd_pulse, largest_time, 
				start_time_btw, end_time_btw, doBTW)

		this_npe = 0
		if excludeFADC:
			this_npe += this_atwd
		elif excludeATWD:
			# print('Omkey {}, portia fadc {}'.format(omkey, this_fadc))
			this_npe += this_fadc
		elif this_fadc > this_atwd:
			this_npe += this_fadc
		elif this_atwd > this_fadc:
			this_npe += this_atwd

		# add the contribution of this dom to the total
		omkey_npe_dict[omkey] = this_npe
		best_npe += this_npe
	
	return best_npe, omkey_npe_dict

def get_digitizer_baselines_dict(calibration, splitted_dom_map, best_pulse_map, atwd_pulse_map, fadc_pulse_map):
	'''
	Get dict of digitizer baselines
	'''
	dict_digitizer_baselines = {}

	for omkey, launch_times in splitted_dom_map:

		this_omkey_baselines = {}
		this_omkey_baselines['beacon'] = []
		# first, the atwd
		# order will be [atwd0 ch0, atwd0 ch1, atwd0 ch2, atwd1 ch0, atwd1 ch1, atwd2 ch2]
		calib = calibration.dom_cal[omkey]
		for atwd in range(2): # first, the atwds
			for channel in range(3):
				beacon_baseline = calib.atwd_beacon_baseline[(atwd,channel)]
				this_omkey_baselines['beacon'].append(beacon_baseline)
		this_omkey_baselines['beacon'].append(calib.fadc_beacon_baseline) # and the fadc


		# then, portia
		# order will be [best, atwd, fadc]
		this_omkey_baselines['portia'] = [-1e15, -1e15, -1e15] # some enormous and unphysical number so I can tag it later
		if omkey in best_pulse_map:
			this_omkey_baselines['portia'][0] = best_pulse_map[omkey].GetBaseLine()
		if omkey in atwd_pulse_map:
			this_omkey_baselines['portia'][1] = atwd_pulse_map[omkey].GetBaseLine()
		if omkey in fadc_pulse_map:
			this_omkey_baselines['portia'][2] = fadc_pulse_map[omkey].GetBaseLine()


		# store the dict
		dict_digitizer_baselines[omkey] = this_omkey_baselines

	return dict_digitizer_baselines

	# # loop over omkeys in the *splitted dom map* (that already passed portia cleaning)
	# for omkey, launches in launch_map:
	# 	if omkey in launch_map:
			
	# 		# loop over launches
	# 		launches = launch_map[omkey]
	# 		for ilaunch, launch in enumerate(launches):
	# 			if launch.lc_bit: # if the lc bit latched

	# 				# get the om calibration
	# 				calib = calibration.dom_cal[omkey]
	# 				atwd_id = launch.which_atwd
	# 				beacon_00 = calib.atwd_beacon_baseline[(0,0)]*I3Units.mV

	# 				print(beacon_00)

	# 				# atwd_id = launch.which_atwd
	# 				# for atwd_chan in range(3):
	# 				# 	baseline = calibration.get_atwd_beacon_baseline(atwd_id, atwd_chan)
	# 				# 	print('Baseline is {}'.format(baseline))						
	# 				# print("Omkey {}, Launch {}, awtd is {}".format(omkey, ilaunch, atwd_id))
	# 				# baseline = calibration.cal.GetATWDBeaconBaseline(atwd_id)


	# 				# trace = launc.raw_atwd[atwd_chan]
	# 				# bins = launch.raw_atwd[atwd_chan]
	# 				# for bin in bins:
	# 				# 	channel = bin.channel
	# 				# 	print('Bin {}, channel {}'.format(bin, channel))


	# 				# # loop over atwd channels (there are 3, 0->2)
	# 				# for atwd_chan in range(3):

	# 				# 	# how many bins did that channel record?
	# 				# 	nbins = len(launch.raw_atwd[atwd_chan])
	# 				# 	if nbins > 0:


def find_high_qe_doms(calibration):

	high_qe_doms = []

	for omkey in calibration.dom_cal.keys():
		dom_cal = calibration.dom_cal[omkey]
		if dom_cal.relative_dom_eff > 1.1:
			high_qe_doms.append(omkey)

	return high_qe_doms

def LoopEHEPulses(frame, doBTW=True, excludeHighQE=False, excludeFADC=False, excludeATWD=False, writeBaselines=False):
	if not frame['I3EventHeader'].sub_event_stream == 'InIceSplit':
		return False
	
	high_qe_doms = []
	if excludeHighQE:
		if not 'I3Calibration' in frame or not 'I3DetectorStatus' in frame:
			icetray.logging.log_fatal('I3Calibration or I3Detector status not in frame, but excludeHighQE=True')

		calibration = frame['I3Calibration']
		high_qe_doms = find_high_qe_doms(calibration)

	# the way it is *actually* calculated in Portia 
	# https://code.icecube.wisc.edu/projects/icecube/browser/IceCube/meta-projects/combo/trunk/portia/private/portia/I3Portia.cxx#L604
	# using the realtime and L2 scripts has makeBestPulse enabled
	# so starting from an L1->L2 file, we actually want to compare EHEFADCPortiaPulses
	# to EHEBestPortiaPulse, where the latter is the "stitched together" FADC
	# and ATWD pulses; if you try and compare to EHEATWDPortiaPulse
	# you do not get the right answer, because the realtime and L2 filters
	# use makeBestPulse=True
	splitted_dom_map = frame.Get('splittedDOMMap')
	fadc_pulse_map = frame.Get('EHEFADCPortiaPulse')
	atwd_pulse_map = frame.Get('EHEATWDPortiaPulse')
	best_pulse_map = frame.Get('EHEBestPortiaPulse')

	which_atwd_pulses = best_pulse_map
	if excludeFADC:
		which_atwd_pulses=atwd_pulse_map

	best_npe, portia_omkey_npe_dict = get_portia_omkey_npe_dict(splitted_dom_map, 
		fadc_pulse_map, which_atwd_pulses, high_qe_doms=high_qe_doms, 
		excludeFADC=excludeFADC, excludeATWD=excludeATWD)


	if writeBaselines:

		if not 'I3Calibration' in frame or not 'I3DetectorStatus' in frame:
			icetray.logging.log_fatal('I3Calibration or I3Detector status not in frame')

		if not 'HLCOfflineCleanInIceRawDataWODC' in frame:
			icetray.logging.log_fatal('HLCOfflineCleanInIceRawDataWODC ism not in the frame')

		calibration = frame['I3Calibration']
		launch_series_map = frame['HLCOfflineCleanInIceRawDataWODC']
		dict_omkey_baselines = get_digitizer_baselines_dict(calibration=calibration, 
			splitted_dom_map = splitted_dom_map,
			best_pulse_map = best_pulse_map, 
			atwd_pulse_map = atwd_pulse_map, 
			fadc_pulse_map = fadc_pulse_map)

		# portia_omkey_baseline_fadc_dict = get_portia_omkey_baseline_dict(splitted_dom_map, 
		# 	fadc_pulse_map)

		# portia_omkey_baseline_best_dict = get_portia_omkey_baseline_dict(splitted_dom_map, 
		# 	best_pulse_map)

		# write to disk
		beacon_atwd_0_0 = []
		beacon_atwd_0_1 = []
		beacon_atwd_0_2 = []
		beacon_atwd_1_0 = []
		beacon_atwd_1_1 = []
		beacon_atwd_1_2 = []
		beacon_fadc = []
		portia_best = []
		portia_atwd = []
		portia_fadc = []
		strings = []
		doms = []
		for omkey, container in dict_omkey_baselines.items():
				strings.append(omkey.string)
				doms.append(omkey.om)
				# beacon
				beacon_atwd_0_0.append(container['beacon'][0])
				beacon_atwd_0_1.append(container['beacon'][1])
				beacon_atwd_0_2.append(container['beacon'][2])
				beacon_atwd_1_0.append(container['beacon'][3])
				beacon_atwd_1_1.append(container['beacon'][4])
				beacon_atwd_1_2.append(container['beacon'][5])
				beacon_fadc.append(container['beacon'][6])
				# portia
				portia_best.append(container['portia'][0])
				portia_atwd.append(container['portia'][1])
				portia_fadc.append(container['portia'][2])
		beacon_atwd_0_0 = np.asarray(beacon_atwd_0_0)
		beacon_atwd_0_1 = np.asarray(beacon_atwd_0_1)
		beacon_atwd_0_2 = np.asarray(beacon_atwd_0_2)
		beacon_atwd_1_0 = np.asarray(beacon_atwd_1_0)
		beacon_atwd_1_1 = np.asarray(beacon_atwd_1_1)
		beacon_atwd_1_2 = np.asarray(beacon_atwd_1_2)
		beacon_fadc = np.asarray(beacon_fadc)
		portia_best = np.asarray(portia_best)
		portia_atwd = np.asarray(portia_atwd)
		portia_fadc = np.asarray(portia_fadc)


		file_out = h5py.File('compare_baselines.hdf5', 'w')
		data_overlap = file_out.create_group('baselines')
		data_overlap.create_dataset('beacon_atwd_0_0', data=beacon_atwd_0_0)
		data_overlap.create_dataset('beacon_atwd_0_1', data=beacon_atwd_0_1)
		data_overlap.create_dataset('beacon_atwd_0_2', data=beacon_atwd_0_2)
		data_overlap.create_dataset('beacon_atwd_1_0', data=beacon_atwd_1_0)
		data_overlap.create_dataset('beacon_atwd_1_1', data=beacon_atwd_1_1)
		data_overlap.create_dataset('beacon_atwd_1_2', data=beacon_atwd_1_2)
		data_overlap.create_dataset('beacon_fadc', data=beacon_fadc)
		data_overlap.create_dataset('portia_best', data=portia_best)
		data_overlap.create_dataset('portia_atwd', data=portia_atwd)
		data_overlap.create_dataset('portia_fadc', data=portia_fadc)
		data_overlap.create_dataset('strings', data=strings)
		data_overlap.create_dataset('doms', data=doms)
		file_out.close()


def get_hit_map(frame,pulse_mask_name):
	if type(frame[pulse_mask_name]) == dataclasses.I3RecoPulseSeriesMap:
		hit_map = frame.Get(pulse_mask_name)
	elif type(frame[pulse_mask_name]) == dataclasses.I3RecoPulseSeriesMapMask:
		hit_map = frame[pulse_mask_name].apply(frame)
	return hit_map

def get_casaulqtot_omkey_npe_dict(calibration, status, vertex_time, 
	causality_window, hit_map, max_q_per_dom,
	exclude_high_qe_doms = True):
	
	causal_qtot=0.;
	omkey_npe_dict = {}; # for all the DOMs
	omkey_npe_dict_noDC = {}; # for DOMs that are not in DC, high QE, or saturated
	
	for omkey, pulses in hit_map.items():
		if not omkey in calibration.dom_cal:
			print("omkey {} not in dom_cal".format(omkey))
			continue
		if not omkey in status.dom_status.keys():
			print("omkey {} not in status".format(omkey))
			continue
		dom_cal = calibration.dom_cal[omkey]

		charge_this_dom = 0.

		# loop pulses
		for p in pulses:
			# print("omkey {}, pulse charge {}".format(omkey, p.charge))
			# print(p.flags)
			if (p.time >= vertex_time and p.time < vertex_time+causality_window):
				charge_this_dom+=p.charge
		omkey_npe_dict[omkey] = charge_this_dom

		# now, to actual compute causal_qtot, we exclude deep core
		# and, exclude any dom with high efficiency
		# and, exclude any DOM which contributes more than max_q_per_dom

		string = omkey.string
		if string in [79, 80, 81, 82, 83, 84, 85, 86]:
			continue
		if dom_cal.relative_dom_eff > 1.1 and exclude_high_qe_doms:
			continue				
		if charge_this_dom < max_q_per_dom:
			omkey_npe_dict_noDC[omkey] = charge_this_dom
			causal_qtot+=charge_this_dom

	return causal_qtot, omkey_npe_dict, omkey_npe_dict_noDC

def LoopHESEPulses(frame, pulses):
	if not frame['I3EventHeader'].sub_event_stream == 'InIceSplit':
		return False

	if not 'I3Calibration' in frame or not 'I3DetectorStatus' in frame:
		icetray.logging.log_fatal('I3Calibration or I3Detector status not in frame')

	calibration = frame['I3Calibration']
	status = frame['I3DetectorStatus']

	vertex_time = frame.Get('HESE_VHESelfVetoVertexTime').value #I3Double is funky
	causality_window = 5000. * I3Units.ns
	hit_map = get_hit_map(frame, pulses)
	causal_qtot, causalqtot_omkey_npe_dict, causalqtot_omkey_npe_dict_noDC = get_casaulqtot_omkey_npe_dict(calibration, 
		status, vertex_time, causality_window, hit_map, Inf)
	
	print("Causal Qtot is {:.2f}".format(causal_qtot))

def Compare_HESE_EHE(frame, hese_pulses, table_name=None, exclude_fadc=False, exclude_atwd=False):

	if not frame['I3EventHeader'].sub_event_stream == 'InIceSplit':
		return False

	# first, HESE
	if not 'I3Calibration' in frame or not 'I3DetectorStatus' in frame:
		icetray.logging.log_fatal('I3Calibration or I3Detector status not in frame')

	calibration = frame['I3Calibration']
	status = frame['I3DetectorStatus']
	
	vertex_time = frame.Get('HESE_VHESelfVetoVertexTime').value #I3Double is funky
	causality_window = 5000. * I3Units.ns
	hit_map = get_hit_map(frame, hese_pulses)
	causal_qtot, causalqtot_omkey_npe_dict,  causalqtot_omkey_npe_dict_noDC= get_casaulqtot_omkey_npe_dict(calibration, 
		status, vertex_time, causality_window, hit_map, Inf)

	# now, EHE
	splitted_dom_map = frame.Get('splittedDOMMap')
	fadc_pulse_map = frame.Get('EHEFADCPortiaPulse')
	atwd_pulse_map = frame.Get('EHEATWDPortiaPulse')
	best_pulse_map = frame.Get('EHEBestPortiaPulse')

	which_atwd_pulses = best_pulse_map
	if exclude_fadc:
		which_atwd_pulses=atwd_pulse_map

	best_npe, portia_omkey_npe_dict = get_portia_omkey_npe_dict(splitted_dom_map, 
		fadc_pulse_map, which_atwd_pulses, excludeFADC=exclude_fadc, excludeATWD=exclude_atwd)

	# check the overlap
	missing_in_hese = []
	missing_in_ehe = []
	num_hese_chans=0
	for omkey, portia_npe in portia_omkey_npe_dict.items():
		if omkey not in causalqtot_omkey_npe_dict_noDC:
			missing_in_hese.append(portia_npe)

	for omkey, hese_npe in causalqtot_omkey_npe_dict_noDC.items():
		if omkey not in portia_omkey_npe_dict:
			missing_in_ehe.append(hese_npe)
	missing_in_hese = np.asarray(missing_in_hese)
	missing_in_ehe = np.asarray(missing_in_ehe)

	# save the overlapping data to an hdf5 file
	hese_overlap = []
	ehe_overlap = []
	strings = []
	doms = []
	for omkey, portia_npe in portia_omkey_npe_dict.items():
		if omkey in causalqtot_omkey_npe_dict_noDC:
			strings.append(omkey.string)
			doms.append(omkey.om)
			hese_npe = causalqtot_omkey_npe_dict_noDC[omkey]
			hese_overlap.append(hese_npe)
			ehe_overlap.append(portia_npe)

	hese_overlap = np.asarray(hese_overlap)
	ehe_overlap = np.asarray(ehe_overlap)
	strings = np.asarray(strings)
	doms = np.asarray(doms)

	if table_name is None:
		table_name='comparison_overlap_{}.hdf5'.format(hese_pulses)
	file_out = h5py.File(table_name, 'w')
	data_overlap = file_out.create_group('data_overlap')
	data_overlap.create_dataset('hese_overlap', data=hese_overlap)
	data_overlap.create_dataset('ehe_overlap', data=ehe_overlap)
	data_overlap.create_dataset('strings', data=strings)
	data_overlap.create_dataset('doms', data=doms)

	data_exclusion = file_out.create_group('data_exclusion')
	data_exclusion.create_dataset('missing_in_hese', data=missing_in_hese)
	data_exclusion.create_dataset('missing_in_ehe',data=missing_in_ehe)
	file_out.close()

@icetray.traysegment
def HeseFilter(tray, name, pulses='RecoPulses', If = lambda f: True):

	from icecube import VHESelfVeto
	# from icecube import clast
	# from icecube import linefit, lilliput

	TriggerEvalList = ['InIceSMTTriggered'] # work on SMT8 triggers
	def If_with_triggers(frame):
		if not If(frame):
			return False
		else:
			# very hacky, the I3 file from pole seems to be missing 
			# from the InIceSMTTriggered bool in the P-frame
			# so I'm shoe-horning in this "true"
			return True
		# for trigger in TriggerEvalList:
		# 	if frame[trigger].value:
		# 		print("Frame value ")
		# 		return True
		# return False

	tray.AddModule('HomogenizedQTot', name+'_qtot_total_brian',
		Pulses=pulses,
		Output='HESE_CausalQTot_Redo',
		If = If_with_triggers)

################################################################
################################################################
################################################################
################################################################

#							GRAVEYARD						   #

################################################################
################################################################
################################################################
################################################################


# essentially create a I3Map<OMKey, I3PulseSeries> from a I3PulseSeriesMapMask
# the former is needed for reconstructors like the HomogenizedQtot of HESE
# pulse_series_name is the name of the map, final_name is the name you want it to hav in the end
def PutPulsesFromPulseMaskIntoFrame(frame, pulse_mask_name, final_name):
	if not frame['I3EventHeader'].sub_event_stream == 'InIceSplit':
		return False
	else:
		if type(frame[pulse_mask_name]) == dataclasses.I3RecoPulseSeriesMap:
			hit_map = frame.Get(pulse_mask_name)
		elif type(frame[pulse_mask_name]) == dataclasses.I3RecoPulseSeriesMapMask:
			hit_map = frame[pulse_mask_name].apply(frame)
		frame.Put(final_name, hit_map)
		return True	
