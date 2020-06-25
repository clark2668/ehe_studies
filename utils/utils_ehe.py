# IceCube import
from icecube import icetray
from icecube.recclasses import I3PortiaEvent
import numpy as np

def has_ehe_objects(frame):
	"""
	Returns true if a P-frame has EHE L2 objects.

	This function will return true only if the p-frame that has been passed
	contains the EHE L2 objects, namely some Ophelia and Portia reconstructions.
	These are generated at the L2 level at pole.
	(see https://code.icecube.wisc.edu/projects/icecube/browser/IceCube/projects/filterscripts/releases/V19-07-00/python/offlineL2/level2_Reconstruction_EHE.py)

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

def get_ophelia_fitqual(frame):
	"""
	Return EHE L2 Ophelia Improved Linefit Fit Quality

	This function will return the fit quality (chi^2/ndf) of the ophelia
	improved linefit

	The EHE L2 filter running online at pole saves the final result
	as an EHEOpheliaSRT_ImpLF frame, which we can get the fit_quality of

	Parameters
	----------
	frame: I3 P-frame containing EHE objects
		An I3 P-frame containing EHE objects

	"""

	fitqual=-10
	if(frame.Has("EHEOpheliaSRT_ImpLF")):
		fitqual = frame.Get("EHEOpheliaSRT_ImpLF").fit_quality
	return fitqual


def get_ophelia_zenith(frame):
	"""
	Return EHE L2 Ophelia zenith

	This function will return the zenith of the event
	as computed by the Portia EHE filter.

	The EHE L2 filter running online at pole saves the final result
	as an I3Particle with the name "EHEOpheliaParticleSRT_ImpFL"
	(see https://code.icecube.wisc.edu/projects/icecube/browser/IceCube/projects/filterscripts/releases/V19-07-00/python/offlineL2/level2_Reconstruction_EHE.py#L76)
	This is, as best I can tell, a pulse series as reconstructed by Portia
	then being passed to improved LineFit

	Because this is an I3Particle, we get it's properties by calling dir.zenith
	(see http://software.icecube.wisc.edu/documentation/projects/dataclasses/particle.html?highlight=i3particle)

	Parameters
	----------
	frame: I3 P-frame containing EHE objects
		An I3 P-frame containing EHE objects

	"""

	zenith=-10
	if(frame.Has("EHEOpheliaParticleSRT_ImpLF")):
		zenith = frame.Get("EHEOpheliaParticleSRT_ImpLF").dir.zenith
	return zenith

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

def get_lognpecut_by_fitqual(fitqual):
	"""
	Return the log10npe cut value by fit quality for the EHE L3 cut

	This function will return log10(npe) cut value for the EHE L3 cut
	based on the fitqual value


	Parameters
	----------
	fitqual: double
		Ophelia fit quality
	"""

	lognpe_cut = 1e30 #ridiculously large number
	if(fitqual>=30):
		if (fitqual<80):
			lognpe_cut = 4.6
		elif (fitqual<120):
			lognpe_cut = 4.6 + 0.015*(fitqual-80)
		else:
			lognpe_cut = 5.2
	return lognpe_cut

def pass_L3(npe, fitqual):
	"""
	Returns true if an event passes the L3 filter

	This function will return true only if the ehe event passes the L3 filter
	This is a fit quality and NPE cut


	Parameters
	----------
	npe: double
		Number of Portia NPE
	fitqual: double
		Ophelia fit quality
	"""

	lognpe_cut = get_lognpecut_by_fitqual(fitqual)
	result=False # start off with failing the cut
	if(np.log10(npe) >= lognpe_cut):
		result = True
	return result

def pass_L2(npe, nchan, fitqual):
	"""
	Returns true if an event passes the L2 filter

	This function will return true only if the ehe event passes the L2 filter
	Which means npe > 25000
				nchan > 100
				fitqual > 30

	Parameters
	----------
	npe: double
		Number of Portia NPE
	nchan: int
		Number of Portia chans
	fitqual: double
		Ophelia fit quality
	"""

	result=False # start off with failing the cut
	if(npe>25000 and nchan>100 and fitqual>30):
		result=True 
	return result








