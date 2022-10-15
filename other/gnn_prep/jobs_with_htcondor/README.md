# How to Get Started

The first thing to do is to activate the "environment" to use. Let's use the python environment we used in the demo. So:

```
source env.sh
```

The main work in creating the list of jobs to do, the so called "dagman file", is done by:

```
python make_dagman.py
```

The file contains several variables which determine the behavior.
In particular, we might pay attention to `input_file_dir` which specifies
the location of the input files.
As well as `output_file_dir`, which specifies where our converted files should go.

(Note that right now, I've hard coded in the 0000000-0000999 part, 
but you might want to make python figure that part out for you!)

If everything is going alright, it should produce a `dagman_run_i3toh5.dag` file.

This file you can then submit to the cluster:

```
condor_submit_dag dagman_run_i3toh5.dag
```

Note, there are a few other tweaks you'll probably need to change!
In particular, the `script_path` in `job.sh` and 
the `log`, `output` and `error` files in `job.sub`.
A few things are defintely going to break, but ask for help, it'll be alright!