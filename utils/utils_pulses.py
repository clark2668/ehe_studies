from icecube import icetray, dataclasses, portia
from icecube.icetray import I3Units
from icecube.phys_services.which_split import which_split
from I3Tray import Inf

import numpy as np

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
		print("Skip adding npe {}".format(npe))
		npe_out = 0
	return npe_out

def get_portia_omkey_npe_dict(splitted_dom_map, fadc_pulse_map, atwd_pulse_map, 
	excluded_doms = [], doBTW=True):
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
		if this_fadc > this_atwd:
			this_npe += this_fadc
		elif this_atwd > this_fadc:
			this_npe += this_atwd

		# add the contribution of this dom to the total
		# omkey_npe_dict[omkey] = this_npe
		best_npe += this_npe
	
	return best_npe, omkey_npe_dict

def CalcPortiaCharge(frame, DOMsToExclude = []):
	if not frame['I3EventHeader'].sub_event_stream == 'InIceSplit':
		return False

	splitted_dom_map = frame.Get('splittedDOMMapSRT')
	fadc_pulse_map = frame.Get('EHEFADCPortiaPulseSRT')
	atwd_pulse_map = frame.Get('EHEATWDPortiaPulseSRT')
	best_pulse_map = frame.Get('EHEBestPortiaPulseSRT')

	best_npe, _ = get_portia_omkey_npe_dict(splitted_dom_map,
		fadc_pulse_map, best_pulse_map, doBTW=True, excluded_doms=DOMsToExclude)

	print('Best NPE {:.2f}'.format(best_npe))