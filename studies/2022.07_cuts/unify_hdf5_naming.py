"""

The point of this script is to unify the naming conventions in the hdf5 file,
so that the hit multiplicity and linefit variables all have the same name for uniform handling.
"""

import h5py

top_dir = "/disk20/users/brian/IceCube/corsika/"
top_dir = "/disk20/users/brian/IceCube/data/"

files_to_fix = [
    # "combo_20787.hdf5",
    # "combo_20787_smaller.hdf5"
    "combo_IC86-I-pass2_L4.hdf5",
    "combo_IC86-II-pass2_L4.hdf5",
    "combo_IC86-III-pass2_L4.hdf5"
]

for f_to_fix in files_to_fix:

    fin = top_dir + f_to_fix
    print("Working on {}".format(fin))

    with h5py.File(fin, 'a') as f:
        top_keys = [key for key in f.keys()]
        
        # fix the linefit name
        if "EHELineFit" in top_keys and 'LineFit' not in top_keys:
            print("    Fixing EHELineFit")
            f["LineFit"] = f["EHELineFit"]
        elif "LineFit_redo" in top_keys and 'LineFit' not in top_keys:
            print("    Fixing LineFit_redo")
            f["LineFit"] = f["LineFit_redo"]
        # elif "LineFit" in top_keys:
            # print("    Fixing LineFit")
            # f["LineFit_unified"] = f["LineFit"]
        
        # now, fix the hit multiplicity label
        if "HitMultiplicityValues" in top_keys and 'CVMultiplicity' not in top_keys:
            print("    Fixing HitMultiplicityValues")
            f["CVMultiplicity"] = f["HitMultiplicityValues"]
        # elif "CVMultiplicity" in top_keys and
            # print("    Fixing CVMultiplicity")
            # f["HitMultiplicityValues_unified"] = f["CVMultiplicity"]            
        
