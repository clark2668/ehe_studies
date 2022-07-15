import numpy as np

import tables
from tables.exceptions import NoSuchNodeError

from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm

def stepped_path(edges, bins):
    """
    Create a stepped path suitable for histogramming
    :param edges: bin edges
    :param bins: bin contents
    """
    if len(edges) != len(bins) + 1:
        raise ValueError("edges must be 1 element longer than bins")

    x = np.zeros((2 * len(edges)))
    y = np.zeros((2 * len(edges)))

    x[0::2], x[1::2] = edges, edges
    y[1:-1:2], y[2::2] = bins, bins
    return x,y

def get_juliet_enu_binning():
    energy_bins = np.logspace(5, 12, 141)
    bin_center = (energy_bins[1:] + energy_bins[:-1]) * 0.5
    return energy_bins, bin_center


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


def make_enu_2d_hist(
        weight_dict,
        n_gen,
        flux,
        var_values=None,
        var_bins=None,
        prop_matrix=None,
        selection_mask=slice(None),
        livetime=1):
    # Calculate the generation solid angle
    # (only required if effective are will be binned in solid angle)
    cos_th_min = weight_dict.col('MinimumCosTheta')[selection_mask]
    cos_th_max = weight_dict.col('MaximumCosTheta')[selection_mask]
    phi_min = weight_dict.col('MinimumPhi')[selection_mask]
    phi_max = weight_dict.col('MaximumPhi')[selection_mask]
    omega_gen = (phi_max - phi_min) * (cos_th_max - cos_th_min)

    if var_values is None:
        if var_bins is None:
            var_values = np.zeros_like(energies)
        else:
            raise ValueError(
                'When providing var bins, you also have to provide ',
                'var values.')
    if var_bins is None:
        var_bins = np.array([-1., 1.])
    delta_phi = np.unique(phi_max - phi_min)[0]
    # int_delta_omega = delta_phi * np.diff(coszen_bins)
    int_delta_omega = omega_gen
    var_values = var_values[selection_mask]

    # Interaction probability
    int_w = weight_dict.col('InteractionWeightW')[selection_mask]

    # Injection area
    injection_r = weight_dict.col('InjectionSurfaceR')[selection_mask]
    injection_area = np.pi * injection_r**2

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
    else:
        e_max = np.power(10, log_e_max)
        e_min = np.power(10, log_e_min)
        int_gen_e = (1 / (1 - gamma) *
            (e_max**(-gamma + 1) -
             e_min**(-gamma + 1)))

    energy_bins, bin_center = get_juliet_enu_binning()
    flux_values = flux(bin_center) * bin_center
    sqm2sqcm = 10000

    if prop_matrix is None:
        # Make a histogram of interaction probabilities
        mask = int_w <= 1
        weights = injection_area * int_gen_e * int_w * omega_gen * livetime

        hist, _ = np.histogram(var_values[mask],
                               bins=[var_bins, energy_bins],
                               weights=weights[mask])
    else:
        n_var_bins = len(var_bins) - 1
        # Get the propagation matrices
        if not isinstance(prop_matrix, list):
            prop_matrices = prop_matrix.col('item').reshape(-1, 140)
        else:
            prop_matrices = prop_matrix[0].col('item').reshape(-1, 140)
            for i in range(len(prop_matrix)-1):
                prop_matrices += prop_matrix[i+1].col('item').reshape(-1, 140)

        prop_matrices = prop_matrices[selection_mask]

        mask = int_w <= 1
        # Make a histogram of interaction probabilities
        weights = (prop_matrices * injection_area[:, np.newaxis] *
                   int_gen_e[:, np.newaxis] * int_w[:, np.newaxis] *
                   omega_gen[:, np.newaxis])[mask]
        weights *= (livetime * sqm2sqcm)

        # Now also bin in var and energy at detector
        indices_var = np.digitize(var_values[mask], var_bins) - 1
        hist = sub_sum_partition(weights, indices_var)
        # The number of bins in hist will be limited by the
        # number of bins with actual entries.
        # Just fill the remaining ones with zeros
        hist = np.pad(
            hist,
            ((0, (n_var_bins)-hist.shape[0]), (0, 0))
        )
        hist *= flux_values[np.newaxis, :]
    area = hist / n_gen

    return area.squeeze()


def calc_effective_area(
        energies,
        weight_dict,
        n_gen,
        energy_bins,
        var_values=None,
        var_bins=None,
        prop_matrix=None,
        selection_mask=slice(None),
        livetime=1):
    # Calculate the generation solid angle
    # (only required if effective are will be binned in solid angle)
    cos_th_min = weight_dict.col('MinimumCosTheta')[selection_mask]
    cos_th_max = weight_dict.col('MaximumCosTheta')[selection_mask]
    phi_min = weight_dict.col('MinimumPhi')[selection_mask]
    phi_max = weight_dict.col('MaximumPhi')[selection_mask]
    omega_gen = (phi_max - phi_min) * (cos_th_max - cos_th_min)

    if var_values is None:
        if var_bins is None:
            var_values = np.zeros_like(energies)
        else:
            raise ValueError(
                'When providing var bins, you also have to provide ',
                'var values.')
    if var_bins is None:
        var_bins = np.array([-1., 1.])
    delta_phi = np.unique(phi_max - phi_min)[0]
    # int_delta_omega = delta_phi * np.diff(coszen_bins)
    int_delta_omega = omega_gen
    var_values = var_values[selection_mask]
    energies = energies[selection_mask]

    # Interaction probability
    int_w = weight_dict.col('InteractionWeightW')[selection_mask]

    # Injection area
    injection_r = weight_dict.col('InjectionSurfaceR')[selection_mask]
    injection_area = np.pi * injection_r**2

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
        weights = injection_area * int_gen_e * int_w * omega_gen * livetime

        hist, _, _ = np.histogram2d(var_values[mask],
                                    energies[mask],
                                    bins=[var_bins, energy_bins],
                                    weights=weights[mask])
    else:
        n_energy_bins = len(energy_bins) - 1
        n_var_bins = len(var_bins) - 1
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
        weights *= livetime

        # Now also bin in var and energy at detector
        indices_var = np.digitize(var_values[mask], var_bins) - 1
        indices_e = np.digitize(energies[mask], energy_bins) - 1
        combined_indices = indices_var + n_coszen_bins * indices_e
        hist = sub_sum_partition(weights, combined_indices)
        # The number of bins in hist will be limited by the
        # number of bins with actual entries.
        # Just fill the remaining ones with zeros
        hist = np.pad(
            hist,
            ((0, (n_var_bins*n_energy_bins)-hist.shape[0]), (0, 0))
        )
        hist = hist.reshape(n_energy_bins, n_var_bins, -1).swapaxes(0, 1)

    area = hist / n_gen
    area /= int_delta_e[np.newaxis, np.newaxis, :]
    area /= int_delta_omega[:, np.newaxis, np.newaxis]

    return area.squeeze()


