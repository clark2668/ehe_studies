# IceCube import
from icecube import icetray

def get_run_subrun_part(filepath):
	"""
	Returns run id, subrun id, and file part for a full I3 file path.

	Parameters
	----------
	filepath: full path to a I3 data file
		full path to a I3 data file

	"""

	# I3 files are structured with the following pathnames:
	# /data/exp/IceCube/2011/filtered/level2pass2a/1113/Run00118920/Level2pass2_IC86.2011_data_Run00118920_Subrun00000000_00000209.i3.zst

	# first, grab everything after the last slash
	last_part = filepath.rsplit('/',1)[-1]

	# then break up on the '_' for parsing
	split = last_part.split('_')

	# split[3] will contain the Run00118920 bit
	# split[4] will contain the Subrun00000000
	# split[5] will contain the 00000209.i3.zst
	run_piece = split[3]
	subrun_piece = split[4]
	part_piece = split[5]

	# to get the run, grab the *last* 8 entries of the run_piece
	# to get the subrun, grab the *last* 8 entries of the subrun_piece
	# to get the part, grab the *first* 8 entries of the part_piece
	run = int(run_piece[-8:])
	subrun = int(subrun_piece[-8:])
	part = int(part_piece[0:8])

	return run, subrun, part