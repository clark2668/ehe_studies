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

# Running Jobs on the Cluster

Okay, now that we've mastered the steps to generating and plotting a single hdf5 file,
we now need to get more statistics. To do that, let's use the computing cluster.

The basic idea of a computer cluster is that there is a "bank" or "cluster"
of computers we can use. So, if we need to do something 100 times, we can just
ask 100 of the computers in the cluster to each do a job for us.

The tool used to manage submission of jobs to the cluster here at MSU is called SLURM.

So, I wrote a demo script that submits jobs to the cluster in the `run_cluster.pbs` file.

You can submit it by doing `sbatch run_cluster.pbs`. Though you will have to change
the parts that are directed to my scripts (baclark) to instead point at your username!

`run_cluster.pbs` will run `create_hdf5.py` on a bunch of simulation files,
and then move the output to the right place. There are more comments in the file.

