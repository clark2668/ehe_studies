from icecube import icetray

def has_ehe_objects(frame):
	passed_cut = False
	passed_cut = (
		frame.Has("EHEOpheliaParticleSRT_ImpLF")
		and frame.Has("EHEOpheliaSRT_ImpLF") 
		and frame.Has("EHEPortiaEventSummarySRT")
	)
	return passed_cut