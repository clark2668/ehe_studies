#!/usr/bin/env python                                                                                                  

from I3Tray import *
import sys, os, glob
import subprocess
from icecube import icetray, dataio, dataclasses
import numpy as np

class FilterEhe(icetray.I3ConditionalModule):
    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddParameter('filter_mask_name',
                          'Name of the filter mask',
                          'FilterMask')
        self.AddParameter('filter_name',
                          'Name of the filter to select',
                          'HighQFilter_17')

    def Configure(self):
        self.filter_mask_name = self.GetParameter(
            'filter_mask_name')
        self.filter_name = self.GetParameter('filter_name')
    
    def Physics(self, frame):
        if self.filter_name.lower() == 'all':
            filter_name = ['EHEFilter_12', 'EHEFilter_13', 'EHEAlertFilter_15',
                           'EHEAlertFilterHB_15', 'HighQFilter_17']
        else:
            filter_name = self.filter_name
        filter_mask = self.filter_mask_name

        if frame.Has(filter_mask):
            if not isinstance(filter_name, list):
                passed = frame[filter_mask][filter_name].condition_passed
                if passed:
                    self.PushFrame(frame)
                    return
            else:
                for filter_name_i in filter_name:
                    passes = []
                    try:
                        passed = frame[filter_mask][filter_name_i].condition_passed
                        passes.append(passed)
                    except KeyError:
                        continue
                if np.sum(passes) >= 1:
                    self.PushFrame(frame)
                    return
        return




@icetray.traysegment
def HighQFilter(tray, name, infiles, output_i3):
    
    files = []
    
    if infiles is not None:
        if isinstance(infiles, str):
            files.append(infiles)
        else:
            for f in infiles:
                print(f)
                files.append(f)
        
        tray.Add(dataio.I3Reader, FilenameList=files)
        
        # Apply EHE Filter
        tray.AddModule(FilterEhe, 'EHE Filter',
                filter_mask_name='FilterMask',
                filter_name='HighQFilter_17')
        
        def FilterHighCharge(frame, cut_value):
            charge = frame['Homogenized_QTot']
            if charge >= cut_value:
                return True
            else:
                return False

        tray.AddModule(FilterHighCharge, 'charge cut',
                   cut_value=(10**3.6),
                #    Streams=[icetray.I3Frame.Physics]
                   )
        
        tray.Add("I3Writer", 
                 filename=output_i3,
                 Streams=[icetray.I3Frame.TrayInfo, icetray.I3Frame.DAQ, icetray.I3Frame.Physics],
                 DropOrphanStreams=[icetray.I3Frame.Calibration, icetray.I3Frame.DAQ]
                 )
        

def make_parser():
    """Make the argument parser"""
    from optparse import OptionParser
    parser = OptionParser()
    
    parser.add_option("-i", "--input", action="store",
    type="string", default="", dest="infile",
    help="Input i3 file(s)  (use comma separated list for multiple files)")

    parser.add_option("-o", "--output", action="store",
        type="string", default="", dest="outfile",
        help="Main i3 output file")

    parser.add_option("--gsiftp", action="store", type="string", default="",
        dest="gsiftp", help="url for gsiftp for file transfer")
    
    return parser

def main(options, stats={}):
    icetray.logging.set_level("WARN")
    
    tray = I3Tray()
    
    tray.AddSegment(HighQFilter,
                    infiles=options['infile'],
                    output_i3=options['outfile']
                    )
    
    tray.Execute()
    
    usagemap = tray.Usage()
    for mod in usagemap:
        print(mod)

### iceprod stuff ###                                                                                                  
iceprod_available = False
try:
    from simprod.modules import ipmodule
except ImportError:
    try:
        from iceprod.modules import ipmodule
    except ImportError:
        print('Module iceprod.modules not found. Will not define IceProd Class')
    else:
        iceprod_available = True
else:
    iceprod_available = True

if iceprod_available:
    TheHighQFilterModule = ipmodule.FromOptionParser(make_parser(),main)
### end iceprod stuff ###

# the business!                                                                                                        
if (__name__=="__main__"):

    from optparse import OptionParser

    parser = make_parser()

    (options,args) = parser.parse_args()

    #------------------------------                                                                                    
    # check paths 1                                                                                                    
    #------------------------------                                                                                    
    if options.infile == "" :
        print("infile is empty. exit now.")
        sys.exit(0)

    #-------------------------------                                                                                   
    # convert options to dictionary                                                                                    
    #-------------------------------                                                                                   
    opts = {}
    for name in parser.defaults:
        value = getattr(options,name)

        if name == 'infile' :
            values = []
            value = value.replace(" ","")
            if ',' in value :
                values = value.split(',') # split into multiple inputs                                                 
            else :
                values.append(value)
            opts[name] = values
        else :
            opts[name] = value

    #------------------------------                                                                                    
    # append gsiftp path if needed                                                                                     
    #------------------------------                                                                                    
    gsiftp = options.gsiftp
    need_to_modify = []
    if gsiftp != "" :
        need_to_modify = ['infile','outfile']

    for key in need_to_modify :
        value = opts[key]
        if value == "" :
            # if value is empty don't add gsiftp path                                                                  
            continue
        elif key == 'infile' :
            mod_infiles = []
            for infilename in value :
                infilepath = gsiftp + infilename
                mod_infiles.append(infilepath)
            opts[key] = mod_infiles
        else :
            newpath = gsiftp + value
            opts[key] = newpath

    # write options (for debug)                                                                                        
    print("===== option list start =====")
    for key, value in opts.items() :
        print("%s is set to : " %(key), value)
    print("===== option list end =====")

    #------------------------------                                                                                    
    # Done. call main function.                                                                                        
    #------------------------------                                                                                    
    main(opts)
