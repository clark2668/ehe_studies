import numpy as np

def get_juliet_weightdict_and_propmatrix(open_file):
    """
    Input: an open hdf5 file, already loaded with tables

    Return the juliet weight dictionary, propagation matrix,
    and number of events per file for that generated set.

    """
    weight_dict = open_file.get_node('/JulietWeightDict')
    prop_matrix = [open_file.get_node(f'/PropagationMatrix{flav}') for flav in ['NuE', 'NuMu', 'NuTau']]
    evts_per_file = weight_dict.col('EventsPerFile')[0]
    return weight_dict, prop_matrix, evts_per_file

def correct_events_per_file(juliet_species, juliet_energy_level):
    n_he = 150
    n_vhe = 20
    nu_scaling = 4
    
    evts_per_file = -1
    
    if juliet_species in ['nue', 'numu', 'nutau']:
        if juliet_energy_level == 'very_high_energy':
            evts_per_file = n_vhe * nu_scaling
        else:
            evts_per_file = n_he * nu_scaling
    else:
        if juliet_energy_level == 'very_high_energy':
            if juliet_species == 'tau':
                evts_per_file = 100
            else:
                evts_per_file = n_vhe
        else:
            evts_per_file = n_he
    
    return evts_per_file

def get_juliet_enu_binning():
    
    '''
    A function to get the energy bins used in juliet simulation
    

    Returns
    -------
    energy_bins, bin_center: numpy.ndarray
        Numpy arrays of the energy bins, and bin centers, 
        used in Juliet simulation.
    
    '''

    energy_bins = np.logspace(5, 12, 141)
    bin_center = (energy_bins[1:] + energy_bins[:-1]) * 0.5
    return energy_bins, bin_center


def get_juliet_prop_matrix(prop_matrices):
    """Returns a single propagation patrices

    Paramters
    ---------
        prop_matrices: hdf5 node
            The node of an hdf5 file corresponding to the multidimensional array
            array of Juliet propagation matrices.
            Can either be one propagation matrix (for one flavor)
            or a list of propagation matrices (for >1 flavor).
            If a list is entered, it returns the sum over the list.
    
    Returns
    -------
        prop_matrix: np.ndarray
            The total propagation matrix, summed across
            flavors. Final size (n_events, 140).
    """
    if not isinstance(prop_matrices, list):
        # if they've passed us a single prop matrix, all we need to do is reshape
        prop_matrix = prop_matrices.col('item').reshape(-1, 140)
    else:
        # otherwise, we need to sum across matrices
        prop_matrix = prop_matrices[0].col('item').reshape(-1, 140)
        for i in range(len(prop_matrices)-1):
            prop_matrix += prop_matrices[i+1].col('item').reshape(-1, 140)
    return prop_matrix


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
    prop_matrix: hdf5 node 
        The node of an hdf5 file corresponding to the multidimensional array
        array of Juliet propagation matrices.
        Can either be one propagation matrix (for one flavor)
        or a list of propagation matrices (for >1 flavor).
        If a list is entered, the list is summed over before being used.
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
    energy_bins, bin_center = get_juliet_enu_binning()
    flux_values = flux(bin_center) * bin_center * np.log(10)
    
    summed_prop_matrix = get_juliet_prop_matrix(prop_matrix) # get the prop matrix as array
    flux_weight = np.sum(flux_values[np.newaxis, :] * summed_prop_matrix,
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
    selection_mask=slice(None),
    livetime=1
):

    energy_bins, bin_center = get_juliet_enu_binning()
    
    # get interaction weights
    int_w = weight_dict.col('InteractionWeightW')
    int_w[int_w > 1] = 0 # eliminate weights > 1 # TODO: Why does this happen?
    int_w = int_w[selection_mask] # apply event selectionm ask

    prop_matrix = get_juliet_prop_matrix(prop_matrix) # get the prop matrix as array
    prop_matrix = prop_matrix[selection_mask] # apply event selection mask

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
    weights = (summed_prop_weight_per_E_bin * flux_values *
        (injection_area * sqm2sqcm * omega_gen * # scaling factors
        livetime / n_gen_per_dec)
    )
    return weights

def sub_sum_partition(a, partition):
    """
    Generalization of np.bincount(partition, a).
    Sums rows of a matrix for each value of array of non-negative ints.

    :param a: array_like
    :param partition: array_like, 1 dimension, nonnegative ints
    :return: matrix of shape ('one larger than the largest value in partition', a.shape[1:]). The i's element is
    the sum of rows j in 'a' s.t. partition[j] == i
    """
    assert partition.shape == (len(a),)
    n = np.prod(a.shape[1:], dtype=int)
    bins = ((np.tile(partition, (n, 1)) * n).T + np.arange(n, dtype=int)).reshape(-1)
    sums = np.bincount(bins, a.reshape(-1))
    if n > 1:
        sums = sums.reshape(-1, *a.shape[1:])
    return sums

def calc_juliet_effective_area(
        energies, weight_dict, n_gen,
        energy_bins,
        coszen_values=None,
        coszen_bins=None,
        prop_matrix=None,
        selection_mask=slice(None)):
    
    '''
    This function returns the (up to three dimensional)
    effective areas, in m^2. 

    Notes
    -----
    If cos(zenith) values and bin  are provide, 
    the returned np.ndarray contains three dimensions:
        first: true cos(zenith)
        second: energy at the detector volume
        third: neutrino energy at the surface
    
    If no cosine zenith values are provided, the two dimensions are:
        second: energy at the detector volume
        third: neutrino energy at the surface
    
    This means that to get the effective area for surface fluxes,
    the user shoudl do area.sum(axis=0), to sum across
    all "energy at detector" bins.
    
    
    Parameters
    ----------
    energies: np.ndarray
        Array of energies of the input particles in GeV
    weight_dict: tables.table.Table
        Tables table containing information from the JulietWeightDict.
    n_gen: int
        Number of generated events for this particle type.
    energy_bin: np.ndarray
        Array of the energy bins into which the effective area is to be binned
        Currently must be 140, and match the juliet energy binning.
    coszen_values: np.ndarray (optional)
        Array of the true cosine zenith values of the particles
    coszen_bins: np.ndarray (optional)
        Array of the true cosine zenith bins in which to calculate the area
    prop_matrix: hdf5 node 
        The node of an hdf5 file corresponding to the multidimensional array
        array of Juliet propagation matrices.
        Can either be one propagation matrix (for one flavor)
        or a list of propagation matrices (for >1 flavor).
        If a list is entered, the list is summed over before being used.
    selection mask: np.mask
        A mask indicating which events should be used in the calculation.
        This can be used to apply cuts to the data.

    Returns
    -------
    area: np.ndarray
        Multi-dimensional array of the effective areas.

    '''
    
    
    # Calculate the generation solid angle
    # (only required if effective are will be binned in solid angle)
    cos_th_min = weight_dict.col('MinimumCosTheta')[selection_mask]
    cos_th_max = weight_dict.col('MaximumCosTheta')[selection_mask]
    phi_min = weight_dict.col('MinimumPhi')[selection_mask]
    phi_max = weight_dict.col('MaximumPhi')[selection_mask]
    omega_gen = (phi_max - phi_min) * (cos_th_max - cos_th_min)

    if coszen_values is None:
        if coszen_bins is None:
            coszen_values = np.zeros_like(energies)
        else:
            raise ValueError(
                'When providing cos(theta) bins, you also have to provide ',
                'cos(theta) values.')
    if coszen_bins is None:
        coszen_bins = np.array([-1., 1.])
    delta_phi = np.unique(phi_max - phi_min)[0]
    int_delta_omega = delta_phi * np.diff(coszen_bins)
    coszen_values = coszen_values[selection_mask]
    energies = energies[selection_mask]

    # Interaction probability
    int_w = weight_dict.col('InteractionWeightW')[selection_mask]

    # Injection area
    r = weight_dict.col('InjectionSurfaceR')[selection_mask]
    injection_area = np.pi * (r**2)

    # Gamma
    gamma = weight_dict.col('PowerLawIndex')[selection_mask]
    if len(np.unique(gamma)) > 1:
        raise ValueError(
            'Different Gammas between different files. '
            'Calc int_delta_e for each event individually!'
            '(Its easier to read this way though.)')
    gamma = np.unique(gamma)[0]

    # Get generation energy range
    log_e_min = weight_dict.col('MinimumEnergyLog')[selection_mask]
    log_e_max = weight_dict.col('MaximumEnergyLog')[selection_mask]

    # Calulate \int_E0^E1 dE \Phi_gen
    # and \int_E^{E+\DeltaE} dE \Phi_gen
    if gamma == 1.:
        int_gen_e = (log_e_max - log_e_min) * np.log(10)
        int_delta_e = np.diff(np.log(energy_bins))
    else:
        e_max = np.power(10, log_e_max)
        e_min = np.power(10, log_e_min)
        int_gen_e = (1 / (1 - gamma) *
            (e_max**(-gamma + 1) -
             e_min**(-gamma + 1)))

        int_delta_e = (1 / (1 - gamma) *
            (energy_bins[1:]**(-gamma + 1) -
             energy_bins[:-1]**(-gamma + 1)))

    if prop_matrix is None:
        # Make a histogram of interaction probabilities
        mask = int_w <= 1
        weights = injection_area * int_gen_e * int_w * omega_gen

        hist, _, _ = np.histogram2d(coszen_values[mask],
                                    energies[mask],
                                    bins=[coszen_bins, energy_bins],
                                    weights=weights[mask])
    else:
        n_energy_bins = len(energy_bins) - 1
        n_coszen_bins = len(coszen_bins) - 1
        # Get the propagation matrices
        if (140 % n_energy_bins) != 0:
            raise ValueError('For now energy bins have to be a divisor '
                             'of the prop_matrix shape!')
        downsample_ratio = 140 // n_energy_bins
        if not isinstance(prop_matrix, list):
            prop_matrices = prop_matrix.col('item').reshape(-1, 140)
        else:
            prop_matrices = prop_matrix[0].col('item').reshape(-1, 140)
            for i in range(len(prop_matrix)-1):
                prop_matrices += prop_matrix[i+1].col('item').reshape(-1, 140)
        prop_matrices = prop_matrices.reshape(
            prop_matrices.shape[0], -1, downsample_ratio
        ).sum(axis=2)

        prop_matrices = prop_matrices[selection_mask]
        bin_centers = (energy_bins[1:] + energy_bins[:-1]) * 0.5

        mask = int_w <= 1
        # Make a histogram of interaction probabilities
        weights = (prop_matrices * injection_area[:, np.newaxis] *
                   int_gen_e[:, np.newaxis] * int_w[:, np.newaxis] *
                   omega_gen[:, np.newaxis])[mask]

        # Now also bin in coszen and energy at detector
        indices_cz = np.digitize(coszen_values[mask], coszen_bins) - 1
        indices_e = np.digitize(energies[mask], energy_bins) - 1
        combined_indices = indices_cz + n_coszen_bins * indices_e
        hist = sub_sum_partition(weights, combined_indices)
        # The number of bins in hist will be limited by the
        # number of bins with actual entries.
        # Just fill the remaining ones with zeros
        hist = np.pad(
            hist,
            ((0, (n_coszen_bins*n_energy_bins)-hist.shape[0]), (0, 0))
        )
        hist = hist.reshape(n_energy_bins, n_coszen_bins, -1).swapaxes(0, 1)

    #int_delta_e = np.tile(int_delta_e, (len(coszen_bins)-1, 1))
    #int_delta_omega = np.tile(int_delta_omega, (len(energy_bins)-1, 1)).T
    area = hist / n_gen
    area /= int_delta_e[np.newaxis, np.newaxis, :]
    area /= int_delta_omega[:, np.newaxis, np.newaxis]

    return area.squeeze()

