from icecube import icetray, dataclasses, portia
from icecube.icetray import I3Units
from icecube.phys_services.which_split import which_split
from I3Tray import Inf

import numpy as np
import operator

def print_waveforms(frame, outputdir, waveform_name, string, dom,
	title_mod=None):

	header = frame['I3EventHeader']
	key = OMKey(string, dom)
	try:
		waveform_series = frame.Get(waveform_name)[key]
	except:
		print('no waveforms found for this dom. skip')
		return

	# make plots
	fig, ax = plt.subplots()
	ax.set_xlabel(r'Time / ns')

	min_time = 0
	for waveform in waveform_series:
		wf_vect = np.array(waveform.waveform) / I3Units.mV
		start_time = waveform.time
		bin_width = waveform.bin_width
		if min_time==0:
			min_time = start_time

		# Skip possible second launches
		if len(wf_vect) == 128:
			if start_time > min_time + 5000:
				continue
		elif len(wf_vect) == 256:
			if start_time > min_time + 20000:
				continue

		time = np.linspace(start_time, start_time + bin_width * len(wf_vect), len(wf_vect))

		ax.plot(time, wf_vect, label=f'{waveform.source}, Ch {waveform.channel}')
		if len(wf_vect) == 128:
			ax.set_ylabel(r'Voltage / mV')
		elif len(wf_vect) == 256:
			ax.set_ylabel('Voltage / mV')

	ax.set_title(f'({string}, {dom}), Ev {header.event_id}.{header.sub_event_id}, {waveform_name} ({title_mod})')

	ax.legend(loc='best')
	if len(wf_vect) == 128:
		ax.set_ylim([-50,7000])
	filename = os.path.join(
		outputdir,
		f'run{header.run_id}_evt{header.event_id}_sevt_{header.sub_event_id}_str{string}_dom{dom}_{title_mod}.png')
	if len(wf_vect) == 256:
		filename = filename.replace('.png', '_fadc.png')
	print(filename)
	fig.savefig(filename)