import os

from icecube import icetray, dataio, dataclasses
from I3Tray import load

load('libjuliet-interface')


in_1 = os.path.join(
    '/data/user/brianclark/IceCube/ehe/output/sim_juliet/raw',
    'mu/high_energy/detector/1/det_00000001.i3.zst'
)
in_2 = os.path.join(
    '/data/user/brianclark/IceCube/ehe/output/sim_juliet/raw',
    'mu/high_energy/L2/1/Level2_00000001.i3.zst'
)
out_path = os.path.join(
    'test_00000001.i3.zst'
)
i_f1 = dataio.I3File(in_1)
i_f2 = dataio.I3File(in_2)
o_f = dataio.I3File(out_path, 'w')

keys = [
    'EventsPerFile',
    'I3JulietParamsTree',
    'I3JulietPrimaryParams',
    'I3JulietPrimaryParticle',
    'JulietRNGState'
]

while i_f2.more():
    frame = i_f2.pop_frame()
    if frame.Stop == icetray.I3Frame.DAQ:
        while i_f1.more():
            mc_truth_frame = i_f1.pop_daq()
            if mc_truth_frame['I3EventHeader'] == frame['I3EventHeader']:
                break
        for key in keys:
            frame[key] = mc_truth_frame[key]
    o_f.push(frame)

i_f1.close()
i_f2.close()
o_f.close()

