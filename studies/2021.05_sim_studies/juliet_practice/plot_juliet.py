from email import utils
import h5py
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import pandas as pd
import tables
import utils_weights

parser = argparse.ArgumentParser()
parser.add_argument("-mu", type=str,
    dest="mu_file", required=True,
    help = "path to mu file"
)
args = parser.parse_args()

n_gen = 10000 * 150
livetime = 86400 * 365 # 1 year
flux_key = 'JulietPropagationMatrixNeutrinoFlux1'

with tables.open_file(args.mu_file) as f:
    weight_dict = f.get_node('/JulietWeightDict')
    print(weight_dict.col('InIceInjectionRectangleArea'))
    # weights = utils_weights.calc_juliet_weight_from_weight_dict(
    #     weight_dict = weight_dict,
    #     n_gen = n_gen,
    #     livetime = livetime,
    #     flux_key = flux_key
    #     )
    # print(weights)






# f = h5py.File(args.mu_file, "r")
# print(type(f['Homogenized_QTot']))
# print(f['JulietWeightDict'].items())

# mu_file = pd.HDFStore(args.mu_file)
# energies = mu_file.get('PolyplopiaPrimary').get('energy')
# print(type(mu_file.get('JulietWeightDict').get('JulietPropagationMatrixNeutrinoFlux1')))


# fig = plt.figure()
# ax = fig.add_subplot(111)
# bins = np.linspace(4, 10, 60)
# ax.hist(np.log10(energies), bins=bins)
# plt.tight_layout()
# fig.savefig('raw_energy.png')