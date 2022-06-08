import simweights
import nuflux
import tables

H3a_SIBYLL23C = ['H3a_SIBYLL23C', 'H3a_SIBYLL23C_pr', 'H3a_SIBYLL23C_conv',
                'H3a_SIBYLL23C_k', 'H3a_SIBYLL23C_K0', 'H3a_SIBYLL23C_K0L',
                'H3a_SIBYLL23C_K0S', 'H3a_SIBYLL23C_pi', 'H3a_SIBYLL23C_mu'
]

nuflux_models = H3a_SIBYLL23C

# a dictionary of custom fluxes
custom_nu_models = {}

def get_neutrino_flux_model(model):
    flux = None
    if model in nuflux_models:
        flux = nuflux.makeFlux(model)
    elif model in custom_nu_models.keys():
        flux = custom_nu_models[model]
    else:
        raise NotImplementedError('Flux model {} is not implemented'.format(model))
    return flux

cr_flux_models = {'GaisserH4a' : simweights.GaisserH4a(),
                    'GaisserH3a' : simweights.GaisserH3a()
                }

def get_cr_flux_model(model):
    flux = cr_flux_models[model] # will raise KeyError, so should be transparent
    return flux

options = ['nugen', 'corsika']

def get_flux_model(model, simtype):
    if simtype is 'nugen':
        flux_model = get_neutrino_flux_model(model)
    elif simtype is 'corsika':
        flux_model = get_cr_flux_model(model)
    else:
        raise NotImplementedError('Simi type {} is not implemented'.format(simtype))
    return flux_model

def get_weighter(open_hdf_file, simtype, nfiles=1):
    '''
    please first open the file by doing

    open_hdf_file = pandas.HDFStore(file, 'r')

    '''
    if simtype is 'nugen':
        weighter = simweights.NuGenWeighter(open_hdf_file, nfiles=nfiles)
    elif simtype is 'corsika':
        weighter = simweights.CorsikaWeighter(open_hdf_file, nfiles=nfiles)
    else:
        raise NotImplementedError('Simi type {} is not implemented'.format(simtype))
    return weighter


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

    injection_area = weight_dict.col('InIceInjectionRectangleArea')
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

    injection_area = weight_dict.col('InIceInjectionRectangleArea')

    weight = (injection_area * sqm2sqcm * omega_gen *
              livetime * int_w * flux_weight / n_gen_per_dec)
    return weight




