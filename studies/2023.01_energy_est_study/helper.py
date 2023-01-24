import numpy as np

def get_scale(est_type=None):

    scales ={
        'charge': 1,
        'millipedeE': 2,
        'truncatedE': 10
    }
    return scales[est_type]

def zenith_cut(reco_zen = None, reco_e = None, est_type = 'charge'):
    '''
    Input: np.ndarrays of values to base the cut on
    Returns: mask, where "True" indicates an event *PASSES* the cut
    '''
    
    scale = get_scale(est_type)
        
    # base cut
    mask1 = npe >= scale*(10**4.6)

    # zenith dependent npe cut above the horizon
    above_horizon = np.cos(reco_zen) >= 0.06
    npe_thresh = np.power(10, 4.6 + 1.85 * np.sqrt(1 - np.power((1-np.cos(reco_zen))/0.94, 2)))
    npe_thresh *= scale
    # these should NOT pass the cuts
    npe_mask = np.logical_and(
        above_horizon,
        npe < npe_thresh)
    
    total_mask = np.logical_and(
        mask1, ~npe_mask
    )
    return total_mask

def get_lognpecut_by_zenith_9yr(zenith, est_type='charge'):
    """
    Return the log10npe cut value by cos(zenith) EHE L4 cut
    This function will return log10(npe) cut value for the EHE L4 cut
    based on the cos(zenith) value
    Parameters
    ----------
    zenith: double
        Ophelia zenith
    """

    scale = get_scale(est_type)
    
    assert abs(zenith)<1.01*np.pi, \
        'zenith angle is larger than pi; check that zenith is in radians'

    lognpe_cut = 1e30 #ridiculously large number
    coszenith = np.cos(zenith)
    if coszenith < 0.06:
        lognpe_cut = np.log10(scale*np.power(10.,4.6))
    else:
        lognpe_cut = np.log10(scale*np.power(10., (4.6 + 1.85*np.sqrt(1-((coszenith-1.)/0.94)**2))))
    return lognpe_cut