import numpy as np


def calc_juliet_weight_from_weight_dict(
        weight_dict,
        n_gen,
        livetime=1,
        flux_key='JulietPropagationMatrixNeutrinoFlux1'):
    '''
    Calculates event weights for Juliet using precalculated flux weights.

    Parameters
    ----------
    weight_dict : tables.table.Table
        Tables table containing information from the JulietWeightDict.
    n_gen: int
        Number of generated events for this particle type.
    livetime: float
        Livetime in seconds.
    flux_key: str
        Name of the precalculated flux key in the Juliet weight dict.

    Returns
    -------
    weight: np.array
        Array of weights for each event.
    '''

    flux_weight = weight_dict.col(flux_key)
    int_w = weight_dict.col('InteractionWeightW')
    int_w[int_w > 1] = 0

    # Calculate the generation solid angle
    cos_th_min = weight_dict.col('MinimumCosTheta')
    cos_th_max = weight_dict.col('MaximumCosTheta')
    phi_min = weight_dict.col('MinimumPhi')
    phi_max = weight_dict.col('MaximumPhi')
    omega_gen = (phi_max - phi_min) * (cos_th_max - cos_th_min)

    sqm2sqcm = 10000
    log_e_min = weight_dict.col('MinimumEnergyLog')
    log_e_max = weight_dict.col('MaximumEnergyLog')
    n_gen_per_dec = n_gen / (log_e_max - log_e_min)

    # injection_area = weight_dict.col('InIceInjectionRectangleArea')
    r = weight_dict.col('InjectionSurfaceR')
    injection_area = np.pi * (r**2.)

    weight = (injection_area * sqm2sqcm * omega_gen *
              livetime * int_w * flux_weight / n_gen_per_dec)
    return weight


def calc_juliet_weight_from_weight_dict_and_prop_matrix(
        weight_dict,
        prop_matrix,
        flux,
        n_gen,
        livetime=1):
    '''
    Calculates event weights for Juliet using an input flux and
    a propagation matrix.

    Parameters
    ----------
    weight_dict: tables.table.Table
        Tables table containing information from the JulietWeightDict.
    prop_matrix: np.array, shape (n_events, 140)
        Propagation Matrix read out via the PropagationMatrixFiller
        from Juliet.
    flux: function
        Function that takes the energy in GeV as an argument and returns
        the (nu+nubar) flux in units of (GeV^-1 sr^-1 s^-1 cm^-2).
    n_gen: int
        Number of generated events for this particle type.
    livetime: float
        Livetime in seconds.

    Returns
    -------
    weight: np.array
        Array of weights for each event.
    '''
    # For now assuming flavor ratio is 1:1:1
    energy_bins = np.logspace(5, 12, 141)
    bin_center = (energy_bins[1:] + energy_bins[:-1]) * 0.5
    flux_values = flux(bin_center) * bin_center * np.log(10)
    flux_weight = np.sum(flux_values[np.newaxis, :] * prop_matrix,
                         axis=1)

    int_w = weight_dict.col('InteractionWeightW')
    int_w[int_w > 1] = 0

    # Calculate the generation solid angle
    cos_th_min = weight_dict.col('MinimumCosTheta')
    cos_th_max = weight_dict.col('MaximumCosTheta')
    phi_min = weight_dict.col('MinimumPhi')
    phi_max = weight_dict.col('MaximumPhi')
    omega_gen = (phi_max - phi_min) * (cos_th_max - cos_th_min)

    sqm2sqcm = 10000
    log_e_min = weight_dict.col('MinimumEnergyLog')
    log_e_max = weight_dict.col('MaximumEnergyLog')
    n_gen_per_dec = n_gen / (log_e_max - log_e_min)

    # injection_area = weight_dict.col('InIceInjectionRectangleArea')
    r = weight_dict.col('InjectionSurfaceR')
    injection_area = np.pi * (r**2.)

    weight = (injection_area * sqm2sqcm * omega_gen *
              livetime * int_w * flux_weight / n_gen_per_dec)
    return weight


def calc_juliet_flux_weight(
    weight_dict,
    prop_matrix,
    flux,
    n_gen,
    livetime=1
):

    energy_bins = np.logspace(5, 12, 141)
    bin_center = (energy_bins[1:] + energy_bins[:-1]) * 0.5

    # get interaction weights
    int_w = weight_dict.col('InteractionWeightW')
    int_w[int_w > 1] = 0 # eliminate weights > 1 # TODO: Why does this happen?

    # weight the propagation matrix by the interaction probability
    # NxM (N events x M energy bins) times list of length N (N events)
    prop_matrix *= int_w[:, np.newaxis] # 
    
    # sum row-wise to get sum of prop matrices.
    # result is array of length 140
    # describing total weight (summed across events)
    # in that energy bin
    summed_prop_weight_per_E_bin = np.sum(prop_matrix, axis=0)

    # evaluate the fluxes at the M energy bins
    flux_values = flux(bin_center) * bin_center * np.log(10)

    # Calculate the generation solid angle
    cos_th_min = weight_dict.col('MinimumCosTheta')[0]
    cos_th_max = weight_dict.col('MaximumCosTheta')[0]
    phi_min = weight_dict.col('MinimumPhi')[0]
    phi_max = weight_dict.col('MaximumPhi')[0]
    omega_gen = (phi_max - phi_min) * (cos_th_max - cos_th_min)

    sqm2sqcm = 10000
    log_e_min = weight_dict.col('MinimumEnergyLog')[0]
    log_e_max = weight_dict.col('MaximumEnergyLog')[0]
    n_gen_per_dec = n_gen / (log_e_max - log_e_min)

    r = weight_dict.col('InjectionSurfaceR')[0]
    injection_area = np.pi * (r**2.)


    # this a list of length  M
    weight = (summed_prop_weight_per_E_bin * flux_values *
        (injection_area * sqm2sqcm * omega_gen * # scaling factors
        livetime / n_gen_per_dec)
    )
    return energy_bins, weight
