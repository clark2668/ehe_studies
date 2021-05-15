from icecube import icetray, dataclasses, portia, recclasses
from icecube.icetray import I3Units
from icecube.phys_services.which_split import which_split
from I3Tray import Inf

import numpy as np
import operator

magsix_strings = [45, 46, 54, 56, 63, 64]
def is_magsix(omkey):
	'''
	Function to identify if a DOM is on a "magnificent six" string
	'''
	result = False
	if omkey.string in magsix_strings:
		if omkey.om >= 33:
			result = True
	return result

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

def get_portia_pulse_values(portia_pulse, largest_time, window_start, 
	window_end, enforce_window):
	npe = portia_pulse.GetEstimatedNPE()
	t10 = portia_pulse.GetRecoPulse().time
	npe_out = npe
	if enforce_window and not ((t10 - largest_time >= window_start) and (t10 - largest_time <= window_end)):
		# print("Skip adding npe {}".format(npe))
		npe_out = 0
	return npe_out

def get_portia_omkey_npe_dict(splitted_dom_map, fadc_pulse_map, atwd_pulse_map, 
	excluded_doms = [], doBTW=True, excludeFADC=False, excludeATWD=False):
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
	start_time_btw = -4400.0#*I3Units.ns
	end_time_btw = 6400.0#*I3Units.ns
	
	best_npe = 0.
	omkey_npe_dict = {}

	for omkey, launches in splitted_dom_map:

		# skip the high qe doms
		if omkey in excluded_doms:
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
			this_npe += this_fadc
		elif this_fadc > this_atwd:
			this_npe += this_fadc
		elif this_atwd > this_fadc:
			this_npe += this_atwd

		# add the contribution of this dom to the total
		omkey_npe_dict[omkey] = this_npe
		best_npe += this_npe
	
	return best_npe, omkey_npe_dict

def CalcPortiaCharge(frame, DOMsToExclude = [], excludeFADC=False, excludeATWD=False):
	best_npe = -20.
	omkey_dict = {}
	if not frame['I3EventHeader'].sub_event_stream == 'InIceSplit':
		return best_npe, omkey_dict

	splitted_dom_map = frame.Get('splittedDOMMapSRT')
	fadc_pulse_map = frame.Get('EHEFADCPortiaPulseSRT')
	atwd_pulse_map = frame.Get('EHEATWDPortiaPulseSRT')
	best_pulse_map = frame.Get('EHEBestPortiaPulseSRT')

	which_atwd_pulses = best_pulse_map
	if excludeFADC:
		which_atwd_pulses=atwd_pulse_map

	best_npe, omkey_dict = get_portia_omkey_npe_dict(splitted_dom_map,
		fadc_pulse_map, best_pulse_map, doBTW=True, excluded_doms=DOMsToExclude,
		excludeFADC=excludeFADC, excludeATWD=excludeATWD)

	return best_npe, omkey_dict

def CalcPortiaCharge_module(frame, DOMsToExclude=[], excludeFADC=False, excludeATWD=False,
	name='PortiaEventSummarySRT'):
	best_npe, omkey_npe_dict = CalcPortiaCharge(frame, DOMsToExclude=DOMsToExclude,
		excludeFADC=excludeFADC, excludeATWD=excludeATWD)
	dummy_portia_event = recclasses.I3PortiaEvent()
	dummy_portia_event.SetTotalBestNPEbtw(best_npe)
	frame.Put(name, dummy_portia_event)



def CalcPortiaCharge_DeepMagSix(frame, DOMsToExclude = [], excludeFADC=False, excludeATWD=False):
	'''
	Calculate Portia charge for the "magnificent six" strings (45, 46, 54, 56, 63, 64)
	In the "deep" region (OM > 33)
	See e.g. https://wiki.icecube.wisc.edu/index.php/Analysis_of_The_Standard_Candle_Luminosity#String-wise_Deep_NPE_distribution
	'''
	portia, portia_dict = CalcPortiaCharge(frame, DOMsToExclude=DOMsToExclude,
		excludeFADC=excludeFADC, excludeATWD=excludeATWD)

	magsix_charge = 0.
	for omkey, q in portia_dict.items():
		if is_magsix(omkey):
			magsix_charge += q
	return magsix_charge

def CalcPortiaCharge_DeepMagSix_module(frame, DOMsToExclude = [], excludeFADC=False, excludeATWD=False,
	name='PortiaEventSummarySRT_DeepMagSix'):
	magsix_charge = CalcPortiaCharge_DeepMagSix(frame, DOMsToExclude=DOMsToExclude,
		excludeFADC=excludeFADC, excludeATWD=excludeATWD)
	dummy_portia_event = recclasses.I3PortiaEvent()
	dummy_portia_event.SetTotalBestNPEbtw(magsix_charge)
	frame.Put(name, dummy_portia_event)


def get_pulse_map(frame,pulse_mask_name):
	if type(frame[pulse_mask_name]) == dataclasses.I3RecoPulseSeriesMap:
		hit_map = frame.Get(pulse_mask_name)
	elif type(frame[pulse_mask_name]) == dataclasses.I3RecoPulseSeriesMapMask:
		hit_map = frame[pulse_mask_name].apply(frame)
	return hit_map

def get_homogqtot_omkey_npe_dict(calibration, status, vertex_time, 
	causality_window, pulse_map, max_q_per_dom,
	exclude_high_qe_doms = True, do_causal = False):
	
	causal_qtot=0.;
	omkey_npe_dict = {}; # for all the DOMs
	omkey_npe_dict_noDC = {}; # for DOMs that are not in DC, high QE, or saturated
	
	for omkey, pulses in pulse_map.items():
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
			# print(p.flags)
			if do_causal:			
				if (p.time >= vertex_time and p.time < vertex_time+causality_window):
					charge_this_dom+=p.charge
			else:
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

def CalcQTot(frame, pulses='SplitInIcePulses', do_causal=False):
	
	qtot = -10.
	qtot_omkey_npe_dict	= {}
	qtot_omkey_npe_dict_noDC = {}

	if not frame['I3EventHeader'].sub_event_stream == 'InIceSplit':
		return qtot, qtot_omkey_npe_dict, qtot_omkey_npe_dict_noDC

	if not 'I3Calibration' in frame or not 'I3DetectorStatus' in frame:
		icetray.logging.log_fatal('I3Calibration or I3Detector status not in frame')

	calibration = frame['I3Calibration']
	status = frame['I3DetectorStatus']

	pulse_map = get_pulse_map(frame, pulses)

	vertex_time = None
	causality_window = None
	if do_causal:
		vertex_time = frame.Get('HESE_VHESelfVetoVertexTime').value
		causality_window = 5000. * I3Units.ns

	qtot, qtot_omkey_npe_dict, qtot_omkey_npe_dict_noDC = get_homogqtot_omkey_npe_dict(calibration, 
		status, vertex_time, causality_window, pulse_map, Inf, do_causal=do_causal)
	qtot, qtot_omkey_npe_dict, qtot_omkey_npe_dict_noDC = get_homogqtot_omkey_npe_dict(calibration, 
		status, vertex_time, causality_window, pulse_map, qtot/2., do_causal=do_causal)

	return qtot, qtot_omkey_npe_dict, qtot_omkey_npe_dict_noDC

def CalcQTot_module(frame, pulses='SplitInIcePulses', do_causal=False):
	CalcQTot(frame, pulses=pulses, do_causal=do_causal)


def CalcQTot_DeepMagSix(frame, pulses='SplitInIcePulses', do_causal=False):
	'''
	Calculate Qtot for the "magnificent six" strings (45, 46, 54, 56, 63, 64)
	In the "deep" region (OM > 33)
	See e.g. https://wiki.icecube.wisc.edu/index.php/Analysis_of_The_Standard_Candle_Luminosity#String-wise_Deep_NPE_distribution
	'''
	hqtot, hqtot_dict, hqtot_dict_noDC = CalcQTot(frame)

	magsix_charge = 0.
	for omkey, q in hqtot_dict.items():
		if is_magsix(omkey):
			magsix_charge += q
	return magsix_charge

def CalcQTOt_DeepMagSix_module(frame, pulses='SplitInIcePulses', do_causal=False,
	name = 'HomogenizedQTot_DeepMagSix'):
	magsix_charge = CalcQTot_DeepMagSix(frame, pulses=pulses, do_causal=do_causal)
	# frame['HomogenizedQtot_DeepMagSix'] = magsix_charge
	frame.Put(name, dataclasses.I3Double(magsix_charge))


def Compare_Portia_QTot(frame):

	portia, portia_dict = CalcPortiaCharge(frame)
	hqtot, hqtot_dict, hqtot_dict_noDC = CalcQTot(frame)
	print('Porta total {:.2f}, HQtot {:.2f}'.format(portia, hqtot))
	print('------')

	# figure out what DOMs are *shared* between the two calculation methods
	shared_doms = set(portia_dict.keys()).intersection(hqtot_dict_noDC.keys())
	
	# store the differences
	diff = {}
	for omkey in shared_doms:
		diff[omkey] =  portia_dict[omkey] - hqtot_dict[omkey]

	# sort by difference (largest to smallest)
	# this seems to work (even though some people online disagree on if dictionaries 
	# are even sortable at a fundamental level...)
	sorted_diff = dict(sorted(diff.items(), key=lambda item: item[1], reverse=True))

	for omkey, _ in sorted_diff.items():
		print('OMKey {}, Portia {:.2f}, HQtot {:.2f}, Diff {:.2f}'.format(omkey,
			portia_dict[omkey], hqtot_dict[omkey], portia_dict[omkey] - hqtot_dict[omkey]))

	print("\n\n")




