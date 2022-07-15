import numpy as np

def get_juliet_e_binning():
    energy_bins = np.logspace(5, 12, 141)
    bin_center = (energy_bins[1:] + energy_bins[:-1]) * 0.5
    return energy_bins, bin_center

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

    energy_bins, bin_center = get_juliet_e_binning()
    
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

def calculate_weights_var_vs_enu(
    weight_dict,
    prop_matrix,
    flux,
    n_gen,
    livetime,
    var_array,
    var_binning
):

    if len(prop_matrix) != len(var_array):
        raise RuntimeError("Length of prop matrix does not match length of variable array")

    energy_bins, bin_center = get_juliet_e_binning()

    # get interaction weights
    int_w = weight_dict.col('InteractionWeightW')
    int_w[int_w > 1] = 0 # eliminate weights > 1 # TODO: Why does this happen?

    # weight the propagation matrix by the interaction probability
    # NxM (N events x M energy bins) times list of length N (N events)
    prop_matrix *= int_w[:, np.newaxis]

    indices_var = np.digitize(var_array, var_binning)
    

    # we now need to make "slices" through the variable axis

from io import StringIO
class AhlersGZKFlux:
    def __init__(self):
        from scipy import interpolate

        logE, logWeight = np.log10(
            np.loadtxt(
                StringIO(
                    """3.095e5  8.345e-13
                       4.306e5 1.534e-12
                       5.777e5 2.305e-12
                       7.091e5 3.411e-12
                       8.848e5 4.944e-12
                       1.159e6 7.158e-12
                       1.517e6 1.075e-11
                       2.118e6 1.619e-11
                       2.868e6 2.284e-11
                       3.900e6 3.181e-11
                       5.660e6 4.502e-11
                       7.891e6 6.003e-11
                       1.042e7 8.253e-11
                       1.449e7 1.186e-10
                       1.918e7 1.670e-10
                       3.224e7 3.500e-10
                       7.012e7 1.062e-9
                       1.106e8 1.892e-9
                       1.610e8 2.816e-9
                       2.235e8 3.895e-9
                       3.171e8 5.050e-9
                       5.042e8 6.529e-9
                       7.787e8 7.401e-9
                       1.199e9 7.595e-9
                       1.801e9 7.084e-9
                       2.869e9 6.268e-9
                       4.548e9 4.972e-9
                       6.372e9 3.959e-9
                       8.144e9 3.155e-9
                       1.131e10 2.318e-9
                       1.366e10 1.747e-9
                       2.029e10 9.879e-10
                       2.612e10 6.441e-10
                       3.289e10 4.092e-10
                       4.885e10 1.828e-10
                       8.093e10 5.691e-11
                       1.260e11 1.677e-11
                       1.653e11 7.984e-12
                       2.167e11 3.631e-12
                       2.875e11 1.355e-12
                       """""
                )
            )
        ).T

        self._interpolant = interpolate.interp1d(
            logE, logWeight + 8, bounds_error=False, fill_value=-np.inf
        )

    def __call__(self, e_center):
        # Return the flux per flavor
        return 10 ** (self._interpolant(np.log10(e_center)) - 8) / e_center ** 2 / 3.