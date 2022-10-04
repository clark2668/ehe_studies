import h5py
import argparse
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("-f", type=str, nargs='+',
    dest="input_files",
    help="paths to input files (absolute or relative)")

args = parser.parse_args()
files = args.input_files

energies_all = []

for file in files:
    file_in = h5py.File(file, "r")
    data = file_in['data']
    energies = np.asarray(data['energies'])
    print(energies)
    # this is NOT the most efficient way to do this, but anyway...
    for e in energies:
        energies_all.append(e)
    file_in.close()


bins = np.logspace(3,8,10)
fig, ax = plt.subplots(1, 1)
ax.hist(energies_all, bins=bins)
ax.set_xlabel("Energy [GeV]")
ax.set_ylabel("Unweighted Counts")
ax.set_yscale('log')
ax.set_xscale('log')
fig.savefig("energy_dist.png")

