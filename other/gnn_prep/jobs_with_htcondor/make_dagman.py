import numpy as np
import os, glob

# we will have a file, a "dagman" file, which contains our list of jobs
# let's give it this name
dag_file_name = 'dagman_run_i3toh5.dag'


# our dagman file must first contain several "intro" lines
# in particular, we must specify the dagman config
instructions = ""
instructions += 'CONFIG config.dagman\n'
with open(dag_file_name, 'w') as f:
    f.write(instructions)


# then, we want to make up lots of jobs, and in each job, 
# the basic idea is to run the i3 -> hdf5 extraction code, which is kep in job.sh

# this is where the i3 files are on disk
input_file_dir = "/data/user/hieule/iceprod/21941/0000000-0000999"

# this is where we want to put the I3 output files
output_file_dir = "/data/user/brianclark/Gen2_optical/demo_output_files/21941/0000000-0000999/"

# get the list of files to runo ver using glob
file_list = sorted(glob.glob(input_file_dir + "/*"))

job_index = 0

# loop over the files
for file in file_list:

    # take the full path to the file and pull off just the filename, e.g.
    # Sunflower_240m_mDOM_2.2x_MuonGun.021941.000999_baseproc.i3.bz2
    basename = os.path.basename(file) 

    # the output file is the same name as the input file, just with .hdf5 appended
    output_file = basename + ".hdf5"

    # now we write the command
    # we have to set the job number, and relevant variables

    instructions = "" # clear the instructions line
    instructions += f'JOB job_{job_index} job.sub \n'
    instructions += f'VARS job_{job_index} inputfile="{file}" outputfile="{output_file_dir}/{output_file}"\n\n'

    # write it to the list of jobs to run
    with open(dag_file_name, 'a') as f:
        f.write(instructions)

    # increment the job counter
    job_index += 1

