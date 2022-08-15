#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/icetray-start
#METAPROJECT /data/user/brianclark/IceCube/EHE/software/build_icetray

import numpy as np
from tqdm import tqdm
from icecube import icetray, dataio, portia, ophelia

from datetime import datetime

# years = np.arange(2010, 2022)

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-y", type=str, nargs='+',
	dest="years",
	help="which year to run over")

args = parser.parse_args()
years = args.years

for y in years:
    print(f"Working on {y}")

    file = f'/home/brianclark/IceCube/ehe_studies/other/icetray_versions/file_list_{y}.txt'
    data = np.genfromtxt(file, delimiter=",", 
        names=['runid', 'path'], dtype=None)

    outfile = f'/home/brianclark/IceCube/ehe_studies/other/icetray_versions/logs_{y}.txt'

    for r, f in tqdm(zip(data['runid'], data['path']), total=len(data['runid'])):

        f = f.decode('UTF-8')

        file_in = dataio.I3File(f)
        
        versions_involved = {}
        
        frameId = 0
        foundQframe = False
        
        # the tray info frame is one we should encounter almost instantly
        while file_in.more() and frameId < 10:
            try:
                frame = file_in.pop_frame()
            except:
                continue
            if frame.Stop.id == 'I':
                # print(f"Frame {frameId} is at ray info frame")

                info_objs = frame.keys()
                for info in info_objs:
                    the_info = frame[info]
                    the_info_compact = str(the_info.print_compact)
                    if 'I3WaveCalibrator' in the_info_compact:

                        # print("Contains the I3 Wave Calibrator!")
                        svn_url = the_info.svn_url
                        versions_involved[info] = svn_url
            if frame.Stop.id == 'Q' and not foundQframe:
                foundQframe = True
                header = frame.Get("I3EventHeader")
                start_time = header.start_time

            frameId+=1
        
        versions_involved_v2 = {}

        for v, url in versions_involved.items():
            dt = datetime.strptime(v, "%Y-%m-%dT%H:%M:%S.%f")
            versions_involved_v2[dt] = url
        
        from collections import OrderedDict
        versions_involved_v2 = OrderedDict(sorted(versions_involved_v2.items()))
        versions_involved_v2 = list(versions_involved_v2.items())
        newest = versions_involved_v2[-1]
        newest_dt = newest[0]
        newest_url = newest[1]

        with open(outfile, 'a') as fout:
            fout.write(f'{r};{start_time};{newest_dt};{newest_url}\n')

