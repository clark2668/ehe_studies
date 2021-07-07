from glob import glob

from I3Tray import I3Tray
from icecube import hdfwriter, simclasses

# original was 21217 (for Kevin's tutorial)
filedir = "/data/sim/IceCube/2016/filtered/level2/neutrino-generator/21220/0000000-0000999/"
files = sorted(glob(filedir + "*.i3.zst"))

tray = I3Tray()
tray.Add("I3Reader", FileNameList=files)
tray.Add(
    hdfwriter.I3HDFWriter,
    SubEventStreams=["InIceSplit"],
    keys=["PolyplopiaPrimary", "I3MCWeightDict"]
    # keys=['I3MCWeightDict', 'I3EventHeader', 'I3MCTree_preMuonProp', 'PolyplopiaPrimary',
    # 'Homogenized_QTot', 'LineFit', 
    # 'EHEOpheliaParticleSRT_ImpLF', 'EHEPortiaEventSummarySRT'],
    output="Level2_IC86.2016_NuMu.21220.hdf5",
)

tray.Execute()
