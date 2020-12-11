# IceCube import
from icecube import icetray

def get_homogenized_qtot(frame):
	"""
	Return homogonized Qtot for this P-frame.


	This is an I3PODHolder ("plain old data") holder, specifically
	a dataclasses.I3Double. For this reason, in order to get the value out
	as a plain double, we need to call .value

	Parameters
	----------
	frame: I3 P-frame containing a Homogonized_QTot object
		An I3 P-frame containing Homogonized_QTot object

	"""

	homogenized_qtot=-10
	if(frame.Has("Homogenized_QTot")):
		print(frame)
		homogenized_qtot = frame.Get("Homogenized_QTot").value
	return homogenized_qtot

def get_causal_qtot(frame):
	"""
	Return causal Qtot for this P-frame.


	This is an I3PODHolder ("plain old data") holder, specifically
	a dataclasses.I3Double. For this reason, in order to get the value out
	as a plain double, we need to call .value

	Parameters
	----------
	frame: I3 P-frame containing a Calcaul_QTot object
		An I3 P-frame containing Causal_QTot object

	"""

	causal_qtot=-10
	if(frame.Has("Causal_QTot")):
		causal_qtot = frame.Get("Causal_QTot").value
	return causal_qtot