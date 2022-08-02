"""

The point of this script is to unify the naming conventions in the hdf5 file,
so that the hit multiplicity and linefit variables all have the same name for uniform handling.
"""


import h5py

top_dir = "/disk20/users/brian/IceCube/unified_naming/original/"

files_to_fix = [
    # "combo_20787.hdf5",
    # "combo_20787_smaller.hdf5",
    # "combo_21218.hdf5",
    # "combo_21220.hdf5",
    # "combo_IC86-III-pass2_L4.hdf5",
    # "combo_IC86-II-pass2_L4.hdf5",
    # "combo_IC86-I-pass2_L4.hdf5",
    # "Level2_prime_00000001.hd5",
    # "mu_high_energy_l2_prime_1k_merged.hd5",
    # "nue_high_energy_l2_prime_1k_merged.hd5",
    # "numu_high_energy_l2_prime_1k_merged.hd5",
    # "nutau_high_energy_l2_prime_1k_merged.hd5",
    # "tau_high_energy_l2_prime_1k_merged.hd5"
    "combo_v2_21218.hdf5"
]

for f_to_fix in files_to_fix:

    fin = top_dir + f_to_fix
    print("Working on {}".format(fin))

    with h5py.File(fin, 'a') as f:
        top_keys = [key for key in f.keys()]
        
        # fix the linefit name
        if "EHELineFit" in top_keys:
            print("    Fixing EHELineFit")
            f["LineFit_unified"] = f["EHELineFit"]
        elif "LineFit_redo" in top_keys:
            print("    Fixing LineFit_redo")
            f["LineFit_unified"] = f["LineFit_redo"]
        elif "LineFit" in top_keys:
            print("    Fixing LineFit")
            f["LineFit_unified"] = f["LineFit"]
        
        # now, fix the hit multiplicity label
        if "HitMultiplicityValues" in top_keys:
            print("    Fixing HitMultiplicityValues")
            f["HitMultiplicityValues_unified"] = f["HitMultiplicityValues"]
        elif "CVMultiplicity" in top_keys:
            print("    Fixing CVMultiplicity")
            f["HitMultiplicityValues_unified"] = f["CVMultiplicity"]            
        
