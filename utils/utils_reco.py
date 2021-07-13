from icecube import icetray, dataclasses
from icecube import linefit, dipolefit, clast, cscd_llh, fill_ratio, tensor_of_inertia,CascadeVariables
from icecube.icetray import I3Units

icetray.load('double-muon', False) #! This has the pulse map splitter module


'''
NB: This is basically just a reduced version of the OfflineCascadeReco
Because I don't need all of it...
https://github.com/icecube/icetray/blob/main/filterscripts/python/offlineL2/level2_Reconstruction_Cascade.py
'''
@icetray.traysegment
def OfflineCascadeReco( tray, name, If = lambda f: True, suffix = '',
                        SRTPulses = '',
                        Pulses = '',
                        TopoPulses = '',
                        CascadeLineFit = 'CascadeLineFit',
                        CascadeDipoleFit = 'CascadeDipoleFit',
                        CascadeLast = 'CascadeLast',
                        CascadeLlhVertexFit = 'CascadeLlhVertexFit',
                        CascadeLlhVertexFitSplit = 'CascadeLlhVertexFitSplit',
                        BadDOMListName = 'BadDomsList',
                        CascadeFillRatio = 'CascadeFillRatio',
                        CascadeSplitPulses = 'CascadeSplitPulses',
                        CascadeLineFitSplit = 'CascadeLineFitSplit',
                        CascadeToISplit = 'CascadeToISplit',
                        CascadeImprovedLineFit = 'CascadeImprovedLineFit',
                        CascadeContainmentTagging = 'CascadeContainmentTagging',
                        ):
    '''
    :param RawPulses:
        Name of the I3RecoPulseSeriesMap to work on. (not cleaned by cascade hit cleaning)
    :param Pulses:
        Name of the I3RecoPulseSeriesMap to work on. (pre-cleaned by cascade hit cleaning)
    :param CascadeLineFit:
        Name of the output linefit fit
    :param CascadeDipoleFit:
        Name of the output dipolefit fit
    :param CascadeLast:
        Name of the output clast fit
    :param CascadeLlhVertexFit:
        Name of the output CascadeLlh fit
    :param BadDOMListName:
        Name of the Bad DOMs list to use in the FillRatio module
        Should be the list that matches the input pulses (HLC pulses with HLC bad dom list, etc.)
    :param CascadeFillRatio:
        Name of the output FillRation fit
    :param CascadeSplitPulses:
        Basename for the average time split cascade pulses
    :param CascadeLineFitSplit:
        Basename for the linefits based on split pulses
    :param CascadeToISplit:
        Basename for the tensor-of-inertia fits based on split pulses
    :param CascadeImprovedLineFit:
        Name of the output improved linefit result based on cascade pulses
    :param suffix:
        Potential suffix to append to the end of all fits, in case of multiple instances
    :param If:
        Python function or module for conditional execution of all fits
    '''

    tray.AddModule( 'I3LineFit', name + '_CascadeLinefit' + suffix,
                    Name = CascadeLineFit + suffix, # ! Name of fit
                    InputRecoPulses = Pulses,
                    LeadingEdge = 'FLE', # ! Use only first leading edge, for Cascades especially
                    If = If,
                    )


    tray.AddModule( 'I3DipoleFit', name + '_CascadeDipolefit' + suffix,
                    AmpWeightPower = 0,
                    DipoleStep = 0,
                    InputRecoPulses = Pulses,
                    MinHits =  5,
                    Name = CascadeDipoleFit + suffix,
                    If = If,
                    )

    tray.AddModule('I3CLastModule', name + '_CascadeLast' + suffix,
                   Name = CascadeLast + suffix,
                   InputReadout = Pulses,
                   If = If,
                   )

    tray.AddModule( 'I3CscdLlhModule', name + '_CascadeLlh' + suffix,
                    InputType = 'RecoPulse', # ! Use reco pulses
                    RecoSeries = Pulses, # ! Name of input pulse series
                    FirstLE = True, # Default
                    SeedWithOrigin = False, # Default
                    SeedKey = CascadeLast + suffix, # ! Seed fit - CLast reco
                    MinHits = 8, # ! Require 8 hits
                    AmpWeightPower = 0.0, # Default
                    ResultName = CascadeLlhVertexFit + suffix, # ! Name of fit result
                    Minimizer = 'Powell', # ! Set the minimizer to use
                    PDF = 'UPandel', # ! Set the pdf to use
                    ParamT = '1.0, 0.0, 0.0, false',   # ! Setup parameters
                    ParamX = '1.0, 0.0, 0.0, false',   # ! Setup parameters
                    ParamY = '1.0, 0.0, 0.0, false',   # ! Setup parameters
                    ParamZ = '1.0, 0.0, 0.0, false',   # ! Setup parameters
                    If = If,
                    )


    tray.AddModule('I3FillRatioModule', name + '_CascadeFillRatio' + suffix,
                   AmplitudeWeightingPower = 0,
                   BadDOMList = BadDOMListName, #! Get the list of Bad DOMs from the frame
                   RecoPulseName = Pulses,
                   ResultName = CascadeFillRatio + suffix,
                   SphericalRadiusMean = 1.0,
                   SphericalRadiusRMS = 1.0,
                   SphericalRadiusMeanPlusRMS = 1.0,
                   SphericalRadiusNCh = 1.0,
                   VertexName = CascadeLlhVertexFit + suffix,
                   If = If,
                   )