from icecube import icetray

def has_ehe_objects(frame):
	"""
	Returns true if a P-frame has EHE L2 objects.

	This function will return true only if the p-frame that has been passed
	contains the EHE L2 objects, namely the Ophelia and Portia objects

	Parameters
	----------
	frame: I3 P-frame
		An I3 P-frame to be checked

	"""
	passed_cut = False
	passed_cut = (
		frame.Has("EHEOpheliaParticleSRT_ImpLF")
		and frame.Has("EHEOpheliaSRT_ImpLF") 
		and frame.Has("EHEPortiaEventSummarySRT")
	)
	return passed_cut

def get_portia_pulses_and_chans(frame):
	"""
	Return EHE L2 Portia NPE and NChans hit.

	This function will return the number of NPE and number of channels
	hit as computed by the Portia EHE filter.

	Parameters
	----------
	frame: I3 P-frame containing EHE objects
		An I3 P-frame containing EHE objects

	"""
	npe=-10
	nchan = -10
	if(frame.Has("EHEPortiaEventSummarySRT")):
		npe = frame.Get("EHEPortiaEventSummarySRT").GetTotalBestNPE()
		nchan = frame.Get("EHEPortiaEventSummarySRT").GetTotalNch()
	return npe, nchan
