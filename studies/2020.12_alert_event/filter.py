from icecube import icetray, dataclasses, portia
from icecube.icetray import I3Units
from icecube.phys_services.which_split import which_split
from I3Tray import Inf

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
		npe_out = npe
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

def LoopPortiaPulses(frame, doBTW=True):
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

	print("Best npe estimate {}".format(best_npe))
	# return omkey_npe_map

def get_hit_map(frame,pulse_mask_name):
	if type(frame[pulse_mask_name]) == dataclasses.I3RecoPulseSeriesMap:
		hit_map = frame.Get(pulse_mask_name)
	elif type(frame[pulse_mask_name]) == dataclasses.I3RecoPulseSeriesMapMask:
		hit_map = frame[pulse_mask_name].apply(frame)
	return hit_map

def calc_q(calibration, status, vertex_time, causality_window, hit_map, max_q_per_dom):
	hqtot=0.;
	for omkey, pulses in hit_map.items():
		if not omkey in calibration.dom_cal:
			print("omkey {} not in dom_cal".format(omkey))
			continue
		if not omkey in status.dom_status.keys():
			print("omkey {} not in status".format(omkey))
			continue
		dom_cal = calibration.dom_cal[omkey]

		# skip deep core by checking both strings and checking dom efficiency
		string = omkey.string
		if string in [79, 80, 81, 82, 83, 84, 85, 86]:
			continue
		charge_this_dom=0.;
		for p in pulses:
			if (p.time >= vertex_time and p.time < vertex_time+causality_window):
				charge_this_dom+=p.charge
		if dom_cal.relative_dom_eff > 1.1:
			print("{}, eff {}, charge {:.2f}".format(omkey, dom_cal.relative_dom_eff, charge_this_dom))
			continue				
		if charge_this_dom < max_q_per_dom:
			hqtot+=charge_this_dom
	return hqtot

def LoopHESEPulses(frame, pulses):
	if not frame['I3EventHeader'].sub_event_stream == 'InIceSplit':
		return False

	if not 'I3Calibration' in frame or not 'I3DetectorStatus' in frame:
		icetray.logging.log_fatal('I3Calibration or I3Detector status not in frame')

	calibration = frame['I3Calibration']
	status = frame['I3DetectorStatus']

	hqtot = 0.
	vertex_time = frame.Get('HESE_VHESelfVetoVertexTime').value #I3Double is funky
	causality_window = 5000. * I3Units.ns
	hit_map = get_hit_map(frame, pulses)
	hqtot = calc_q(calibration, status, vertex_time, causality_window, hit_map, Inf)
	# hqtot = calc_q(calibration, status, vertex_time, causality_window, hit_map, hqtot/2.)
	print("Causal qtot including high eff DOMs is {:.2f}".format(hqtot))


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
