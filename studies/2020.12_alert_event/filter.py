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

def get_portia_omkey_npe_dict(splitted_dom_map, fadc_pulse_map, atwd_pulse_map, doBTW=True):
	'''
	This is a roughly replication of the MakePortiaEvent function in I3Portia.xx
	We check all of the launched omkeys in the splitted_dom_map
	And get the FADC and ATWD pulse values for that omkey
	And first figure out if they are inside a basetime window (BTW).
	If they are, we figure out which is larger--the FADC or the ATWD--and use that
	as the NPE estimate for that event.
	We will return a dict of the omkey to the NPE estimate (the omkey_npe_dict)
	'''

	# now set up window for cleaning
	largest_time_fadc = find_time_of_largest_pulse(fadc_pulse_map)
	largest_time_best = find_time_of_largest_pulse(atwd_pulse_map)
	largest_time = max(largest_time_fadc, largest_time_best)
	start_time_btw = -4400.0*I3Units.ns
	end_time_btw = 6400.0*I3Units.ns
	
	best_npe = 0.
	omkey_npe_dict = {}

	for omkey, launches in splitted_dom_map:
		
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
		if this_fadc > this_atwd:
			this_npe += this_fadc
		elif this_atwd > this_fadc:
			this_npe += this_atwd

		# add the contribution of this dom to the total
		omkey_npe_dict[omkey] = this_npe
		best_npe += this_npe
	
	return best_npe, omkey_npe_dict

def LoopEHEPulses(frame, doBTW=True):
	if not frame['I3EventHeader'].sub_event_stream == 'InIceSplit':
		return False

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
	best_pulse_map = frame.Get('EHEBestPortiaPulse')

	best_npe, portia_omkey_npe_dict = get_portia_omkey_npe_dict(splitted_dom_map, 
		fadc_pulse_map, best_pulse_map)

	print("Portia NPE is {:.2f}".format(best_npe))
	# return omkey_npe_map

def get_hit_map(frame,pulse_mask_name):
	if type(frame[pulse_mask_name]) == dataclasses.I3RecoPulseSeriesMap:
		hit_map = frame.Get(pulse_mask_name)
	elif type(frame[pulse_mask_name]) == dataclasses.I3RecoPulseSeriesMapMask:
		hit_map = frame[pulse_mask_name].apply(frame)
	return hit_map

def get_casaulqtot_omkey_npe_dict(calibration, status, vertex_time, causality_window, hit_map, max_q_per_dom):
	
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
			if (p.time >= vertex_time and p.time < vertex_time+causality_window):
				charge_this_dom+=p.charge
		omkey_npe_dict[omkey] = charge_this_dom

		# now, to actual compute causal_qtot, we exclude deep core
		# and, exclude any dom with high efficiency
		# and, exclude any DOM which contributes more than max_q_per_dom

		string = omkey.string
		if string in [79, 80, 81, 82, 83, 84, 85, 86]:
			continue
		if dom_cal.relative_dom_eff > 1.1:
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

def Compare_HESE_EHE(frame, hese_pulses):
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
	best_pulse_map = frame.Get('EHEBestPortiaPulse')

	best_npe, portia_omkey_npe_dict = get_portia_omkey_npe_dict(splitted_dom_map, 
		fadc_pulse_map, best_pulse_map)


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

	file_out = h5py.File('comparison_overlap.hdf5', 'w')
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
