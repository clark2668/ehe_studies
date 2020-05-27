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
		homogenized_qtot = frame.Get("Homogenized_QTot").value
	return homogenized_qtot