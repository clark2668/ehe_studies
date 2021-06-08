# How to Get Started

The first thing to do is to activate the "environment" to use. Recall from yesterday:

```
eval `/cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/setup.sh`
$SROOT/metaprojects/combo/V01-01-00/env-shell.sh
```

Then, the first thing we need to do is get all of the data out of I3 format and into a format that's friendlier to plotting, like HDF5. Do to this, we can run the `create_hdf5.py` code. It take two arguments. The file we want to run on, given with the `-i` flag, and the name of the output file we want, give with the `-o` flag.

```
python create_hdf5.py -i /path/to/file.i3 -o output_name
```

and should produce a `output_name.hdf5` file.

Then, we can take that file as an input and read the results back out. This is where you'd imagine making your plots. You can give it the list of files to run over by using the `-f` flag.

```
python plot_variables.py -f output_name.hdf5
```
