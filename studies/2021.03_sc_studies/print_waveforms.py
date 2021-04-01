#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division, print_function

from icecube import dataio, dataclasses, icetray
from I3Tray import *
from icecube.icetray import I3Units
from icecube.icetray import OMKey

import numpy as np
import os

import matplotlib
matplotlib.use('pdf')
from matplotlib import pyplot as plt


WAVEFORM_THRESH = 480000


class plot_waveforms(icetray.I3ConditionalModule):
    def __init__(self, ctx):
        icetray.I3ConditionalModule.__init__(self, ctx)
        self.AddParameter('key',
                          'waveform_key',
                          'CalibratedWaveformsHLCATWD')
        self.AddParameter('outputfolder',
                          'outputfolder',
                          '')
        self.AddParameter('name', 'name', '')
        self.AddParameter('string', 'string', None)
        self.AddParameter('dom', 'dom', None)

    def Configure(self):
        self.key = self.GetParameter('key')
        self.outputfolder = self.GetParameter('outputfolder')
        self.out_name = self.GetParameter('name')
        self.string = self.GetParameter('string')
        self.dom = self.GetParameter('dom')

    def DAQ(self, frame):
        evt_header = frame['I3EventHeader']
        self.run_id = evt_header.run_id
        self.evt_id = evt_header.event_id

    def Physics(self, frame):
        print(self.run_id, self.evt_id)

        wfs = frame[self.key]
        if self.string is None and self.dom is None:
            prev_string = 1
            wf_cnt = 0
            fig, ax = plt.subplots()
            for omkey in wfs.keys():
                if omkey.string != prev_string:
                    prev_string = omkey.string
                    if wf_cnt > 0:
                        ax.legend(loc='best')
                        ax.set_xlabel(r'Time / ns')
                        ax.set_ylabel(r'ATWD Voltage / mV')

                        filename = os.path.join(
                            self.outputfolder,
                            self.out_name + '_evt_{}.pdf')
                        exists = True
                        i = 0
                        while exists:
                            fname = filename.format(i)
                            if os.path.isfile(fname):
                                i += 1
                            else:
                                exists = False

                        print(fname)
                        fig.savefig(fname)
                        wf_cnt = 0
                        fig, ax = plt.subplots()

                waveform_series = wfs[omkey]
                for waveform in waveform_series:
                    wf_vect = np.array(waveform.waveform) / I3Units.mV
                    if np.sum(wf_vect) > WAVEFORM_THRESH:
                        start_time = waveform.time
                        bin_width = waveform.bin_width

                        time = np.linspace(start_time, start_time + bin_width * 128, 128)
                        ax.plot(time, wf_vect, label=omkey)
                        wf_cnt += 1

        elif self.string is not None and self.dom is not None:
            fig, ax = plt.subplots()
            key = OMKey(self.string, self.dom)
            ax.legend(loc='best')
            ax.set_xlabel(r'Time / ns')

            try:
                waveform_series = wfs[key]
            except KeyError:
                print('No waveforms found for this dom. Skipping this event!')
                return
            else:
                min_time = 0
                for waveform in waveform_series:
                    wf_vect = np.array(waveform.waveform) / I3Units.mV
                    start_time = waveform.time
                    bin_width = waveform.bin_width
                    if min_time == 0:
                        min_time = start_time

                    # Skip possible second launches
                    if len(wf_vect) == 128:
                        if start_time > min_time + 5000:
                            continue
                    elif len(wf_vect) == 256:
                        if start_time > min_time + 20000:
                            continue

                    time = np.linspace(start_time, start_time + bin_width * len(wf_vect), len(wf_vect))
                    ax.plot(time, wf_vect, label=key)
                    if len(wf_vect) == 128:
                        ax.set_ylabel(r'ATWD Voltage / mV')
                    elif len(wf_vect) == 256:
                        ax.set_ylabel('FADC Voltage / mV')

                filename = os.path.join(
                    self.outputfolder,
                    self.out_name + f'evt{self.evt_id}_str{self.string}_dom{self.dom}.png')
                if len(wf_vect) == 256:
                    filename = filename.replace('.png', '_fadc.png')
                print(filename)
                fig.savefig(filename)

        else:
            raise NotImplementedError()


def main(inputfile, outputfolder, name, string, dom):
    tray = I3Tray()
    # Read files
    tray.AddModule('I3Reader', 'reader', Filenamelist=[inputfile])

    key = 'CalibratedWaveforms_ATWD'

    tray.AddModule(plot_waveforms, 'plot_wfs',
                   key=key,
                   outputfolder=outputfolder,
                   name=name,
                   string=string,
                   dom=dom)

    tray.AddModule('TrashCan', 'YesWeCan')
    tray.Execute()
    tray.Finish()


if __name__ == '__main__':
    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option('-i', '--infile', type='string', dest='in_file')
    parser.add_option('-o', '--outputdir', type='string', dest='out_dir')
    parser.add_option('-n', '--name', type='string', dest='name')
    parser.add_option('-s', '--string', type=int, dest='string', default=None)
    parser.add_option('-d', '--dom', type=int, dest='dom', default=None)
    (options, args) = parser.parse_args()

    in_file = options.in_file
    out_dir = options.out_dir
    name = options.name
    string = options.string
    dom = options.dom

    main(in_file, out_dir, name, string, dom)

