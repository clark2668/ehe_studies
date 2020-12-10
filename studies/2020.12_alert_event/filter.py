from icecube import icetray, dataclasses
from icecube.icetray import I3Units
from icecube.phys_services.which_split import which_split
from I3Tray import Inf


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
			# continue				
		if charge_this_dom < max_q_per_dom:
			hqtot+=charge_this_dom
	return hqtot

def LoopPulses(frame, pulses):
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
