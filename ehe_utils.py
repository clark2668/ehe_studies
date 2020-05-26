from icecube import icetray

def has_ehe_objects(frame):
	"""
	Return rue if a P-frame has EHE L2 objects.

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
