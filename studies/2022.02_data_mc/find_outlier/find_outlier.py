import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import argparse


# path_to_files = '/data/user/brianclark/IceCube/EHE/l2/20787/*/*'
# path_to_files = '/data/user/brianclark/IceCube/EHE/l2/20787/0000000-0000999/*000220*'
path_to_files = '/data/user/brianclark/IceCube/EHE/l2/20787/0000000-0000999/*'
import glob
file_list = sorted(glob.glob(path_to_files))

for file in file_list:

    cor_file = pd.HDFStore(file)
    # print(cor_file.get('I3EventHeader'))
    try:
        npe = np.asarray(cor_file.get('Homogenized_QTot').get('value'))
        speed = np.asarray(cor_file.get('LineFit').get('speed'))
        evids = np.asarray(cor_file.get('I3EventHeader').get('Event'))
        runidds = np.asarray(cor_file.get('I3EventHeader').get('Run'))
        # # print(npe)
        # # print(cor_file.get('LineFit'))
        # exists = np.asarray(cor_file.get('LineFit').get('exists'))
        # for ev, s, e in zip(evids, speed, exists):
        #     print("Ev {}, Speed {}, Exists {}".format(ev, s, e))

        for n, s, e, r in zip(npe, speed, evids, runidds):
            if np.log10(n) > 4.4 and s < 0.1:
                print("File {}, Run {}, Ev {}, HQtot {}, Speed {}".format(file, r, e, n, s))
        cor_file.close()
    except:
        continue
