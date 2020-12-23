# How to Get Started

The first thing to do is to activate the "environment" to use. I've written a `env.sh` file that contains all the relevant commands. Just do:

```
source env.sh
```

Then, the first thing we need to do is get all of the data out of I3 format and into a format that's friendlier to plotting, like HDF5. Do to this, we can run the `get_from_i3.py` code. It take three arguments. The file we want to run on, given with the `-i` flag, the directory where we want the output writte, given by the `-o` flag, and if we want it to print it's progress "verbosely", by using the `-v` flag.

```
python get_from_i3.py -i /path/to/file.i3 -o /path/to/output/location -v 1
```

If everything is going alright, it should print it's progress like this:

```
Particle     0: NPE=    750.93, Zenith = 1.465
Particle     1: NPE=    394.67, Zenith = 0.000
Particle     2: NPE=    814.01, Zenith = 0.415
Particle     3: NPE=    583.88, Zenith = 0.955
Particle     4: NPE=    547.15, Zenith = 0.628
```

and should produce a `output.hdf5` file.

Then, we can take that file as an input and read the results back out. This is where you'd imagine making your plots. You can give it the list of files to run over by using the `-f` flag.

```
python plot_from_hdf5.py -f output.hdf5
```