# How to Get Started

The first thing to do is to activate the "environment" to use. I've written a `env.sh` file that contains all the relevant commands. Just do:

```
source env.sh
```

Then, the first thing we need to do is get all of the data out of I3 format and into a format that's friendlier to plotting, like HDF5. Do to this, we can run the `get_from_i3.py` code. It take three arguments. The file we want to run on, given with the `-i` flag, the path to the file where we want the output written, given by the `-o` flag.

```
python get_from_i3.py -i /path/to/file.i3 -o /path/to/output/file.hdf5
```

If everything is going alright, it should produce a `file.hdf5` file.

Then, we can take that file as an input and read the results back out. This is where you'd imagine making your plots. You can give it the list of files to run over by using the `-f` flag.

```
python plot_from_hdf5.py -f output.hdf5
```

I have put a sample data file on the HPCC at this location: `/mnt/research/IceCube/Gen2GNN/training_data_from_hieu/21941/`. See if you can start with an I3 file in that location.