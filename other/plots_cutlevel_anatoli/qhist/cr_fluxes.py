import simweights
from simweights._fluxes import CosmicRayFlux
from simweights._pdgcode import PDGCode
import numpy as np
from numpy import exp
from scipy.interpolate import CubicSpline


SEC_PER_YR = 60 * 60 * 24 * 365
# digitized fig.2 of https://pos.sissa.it/395/368/pdf

# Auger protons
auger_p = np.array([
        18.05400, 3.034142181137292e+37,
        18.15659, 2.774046728597764e+37,
        18.24838, 2.501517772832688e+37,
        18.33747, 2.1200558666454891e+37,
        18.42657, 1.712115007694475e+37,
        18.51566, 1.326644286635897e+37,
        18.56425, 1.1243414031991027e+37,
        18.62635, 8.652184035820504e+36,
        18.68035, 6.61240250041143e+36,
        18.75324, 4.4330903924712106e+36,
        18.81533, 2.6071598815378316e+36,
        18.86393, 1.62024555367527e+36,
        18.92603, 7.537815600410613e+35,
        18.95572, 5.0884654064764654e+35

]).reshape(-1, 2)


# Auger A=2-4
auger_he = np.array([
        17.45, 8e35,
        17.47, 9e35,
        17.75, 2.2e36,
        18.04860, 5.8407291303787e+36,
        18.11069, 7.182694997442272e+36,
        18.16199, 8.41685219168439e+36,
        18.23218, 1.049441136383304e+37,
        18.34557, 1.4711761328754588e+37,
        18.44816, 1.8726434238654274e+37,
        18.55076, 2.179331786598341e+37,
        18.64795, 2.271367593021487e+37,
        18.74514, 2.1944079362993673e+37,
        18.85043, 2.034151210812479e+37,
        18.95302, 1.8343113524438305e+37,
        19.05022, 1.6314540647748195e+37,
        19.15011, 1.326644286635897e+37,
        19.25000, 7.965216868160519e+36,
        19.30400, 4.5885624490470074e+36,
        19.37149, 2.2870804266840563e+36,
        19.41199, 1.281694235747813e+36,
        19.44708, 8.30159757086222e+35,
        19.47678, 5.123666318374453e+35

]).reshape(-1, 2)


# Auger A=5-22
auger_n = np.array([
        18.025, 8e35,
        18.05130, 9.079957846749099e+35,
        18.11339, 1.0937602996582911e+36,
        18.18089, 1.3543675101108155e+36,
        18.27268, 1.847000724809125e+36,
        18.35367, 2.518822736058035e+36,
        18.46706, 3.888842544002964e+36,
        18.54806, 5.303357429811899e+36,
        18.62905, 7.332793603807374e+36,
        18.72894, 1.0713715750403207e+37,
        18.80454, 1.4115640814419394e+37,
        18.87203, 1.80919354657282e+37,
        18.98002, 2.484331699177394e+37,
        19.10151, 3.3415817918400445e+37,
        19.25540, 4.1093435911891807e+37,
        19.35799, 4.2534615862528017e+37,
        19.44978, 4.1093435911891807e+37,
        19.54968, 3.506793392946488e+37,
        19.64687, 2.467263698290118e+37,
        19.74676, 1.326644286635897e+37,
        19.80616, 7.91049376954284e+36,
        19.86285, 4.684450806977624e+36,
        19.91955, 2.661642365470281e+36,
        19.96004, 1.7239590519653123e+36,
        20.01674, 9.661174805176381e+35,
        20.07073, 5.642830511827538e+35

]).reshape(-1, 2)


# Auger A=23-38
auger_si = np.array([
        18.4, 2.5e35,
        18.47, 4e35,
        18.53186, 5.1591107426282235e+35,
        18.57775, 6.3008798408642344e+35,
        18.62635, 7.748569930604753e+35,
        18.66955, 9.333829861617642e+35,
        18.72354, 1.1718236949614955e+36,
        18.78024, 1.501919675144877e+36,
        18.82343, 1.82170915934639e+36,
        18.87743, 2.3029019583585864e+36,
        18.92873, 2.85160797328794e+36,
        19.00432, 3.942833004154797e+36,
        19.09611, 5.760750177925207e+36,
        19.17441, 7.856146632756695e+36,
        19.25810, 1.0640109754072502e+37,
        19.34989, 1.4018662713618781e+37,
        19.45518, 1.7478935109837316e+37,
        19.56048, 2.0062969345156018e+37,
        19.67927, 2.0623921999107401e+37,
        19.75486, 1.9652290427820874e+37,
        19.85205, 1.60911404794404e+37,
        19.95194, 1.0640109754072502e+37,
        20.00864, 6.891652222854776e+36,
        20.06263, 4.4026339156742e+36,
        20.10313, 2.871334785404526e+36,
        20.14903, 1.677068872418216e+36,
        20.18413, 7.965216868160518e+35,
        20.20572, 5.0884654064764654e+35

]).reshape(-1, 2)


# Auger A>38
auger_fe = np.array([
        18.78, 1.25e35,
        18.85, 2e35,
        19.01512, 5.194800363804006e+35,
        19.08261, 7.182694997442155e+35,
        19.15281, 9.795305016784596e+35,
        19.22030, 1.3450626552591373e+36,
        19.30940, 2.0910252707092255e+36,
        19.39039, 2.9925947120070834e+36,
        19.45518, 3.997573036899032e+36,
        19.54698, 5.56556144778882e+36,
        19.64687, 7.232383340311357e+36,
        19.74946, 8.244563449874323e+36,
        19.84935, 8.533706941579005e+36,
        19.94384, 7.965216868160519e+36,
        20.04914, 6.891652222854776e+36,
        20.14633, 5.159110742628223e+36,
        20.19762, 4.1093435911891135e+36,
        20.26242, 2.871334785404526e+36,
        20.29752, 2.1494895067305476e+36,
        20.34611, 1.4711761328754588e+36

]).reshape(-1, 2)

auger_components = [auger_p, auger_he, auger_n, auger_si, auger_fe]
auger_comp_names = ['PPlus', 'He4Nucleus', 'N14Nucleus', 'Al27Nucleus', 'Fe56Nucleus']
auger_splines = [CubicSpline(comp[:, 0], comp[:, 1], extrapolate=False) for comp in auger_components]

def auger_flux(energies):
    # trying to create fig4 from https://pos.sissa.it/395/324/pdf
    phi_0 = 8.34e-11 # km^-2 sr-^1 yr^-1 eV^-1
    e_pivots = np.array([
        1e16,    # E_pivot
        2.8e16,  # low energy ankle
        1.58e17, # 2nd knee
        5.0e18,  # ankle
        1.4e19,  # instep
        4.7e19   # suppression
    ])
    gammas = np.array([
        3.09,
        2.85,
        3.283,
        2.54,
        3.03,
        5.3
    ])
    omegas = np.array([
        0.25,
        0.25,
        0.05,
        0.05,
        0.05
    ])
    prod_term = (1 + (energies[:, np.newaxis]/e_pivots[np.newaxis, 1:])**(1/omegas[np.newaxis, :])
                )**(-omegas[np.newaxis, :]*np.diff(gammas)[np.newaxis, :])

    return phi_0 * (energies/e_pivots[0])**(-gammas[0]) * np.prod(prod_term,axis=1)


class GaisserH4aAugerZombie(CosmicRayFlux):
    r"""
    Spectral fits from an internal report\ [#Gaisser1]_ and in Astropart. Phys\ [#Gaisser2]_ by Tom Gaisser.
    In the model H4a, on the other hand, the extra-galactic component
    is assumed to be all protons.
    """
    pdgids = [
        PDGCode.PPlus,
        PDGCode.He4Nucleus,
        PDGCode.N14Nucleus,
        PDGCode.Al27Nucleus,
        PDGCode.Fe56Nucleus,
    ]

    def _flux(idx):
        _gaisser_funcs = [
            lambda E: (0.7860 * E**-2.66 * exp(-E / (4e6 * 1)) +
                       0.0020 * E**-2.4 * exp(-E / (3e7 * 1)) +
                       0.0200 * E**-2.6 * exp(-E / 6e10)),
            lambda E: (0.3550 * E**-2.58 * exp(-E / (4e6 * 2)) +
                       0.002 * E**-2.4 * exp(-E / (3e7 * 2))),
            lambda E: 0.2200 * E**-2.63 * exp(-E / (4e6 * 7)) + 0.00134 * E**-2.4 * exp(-E / (3e7 * 7)),
            lambda E: 0.1430 * E**-2.67 * exp(-E / (4e6 * 13)) + 0.00134 * E**-2.4 * exp(-E / (3e7 * 13)),
            lambda E: 0.2120 * E**-2.63 * exp(-E / (4e6 * 26)) + 0.00134 * E**-2.4 * exp(-E / (3e7 * 26)),
        ]
        _e_breaks = [
            10**(auger_comp[0, 0] - 9) for auger_comp in auger_components
        ]

        def f(E):
            gaisser_vals = _gaisser_funcs[idx](E)
            auger_vals = auger_splines[idx](np.log10(E*1e9))/((E*1e9)**3)*1e9/SEC_PER_YR/1e10
            m = E >= _e_breaks[idx]
            # m1 = np.isfinite(gaisser_vals)
            # m2 = np.isfinite(auger_vals)
            # m = np.logical_and(m1, m2)
            # minimum_overlapped_E = np.min(E[m])
            # E_idx = np.where(E == minimum_overlapped_E)[0][0]
            flux = gaisser_vals
            flux[m] = auger_vals[m]
            return flux
        return f

    _funcs = [_flux(0), _flux(1), _flux(2), _flux(3), _flux(4)]


def main():
    from matplotlib import pyplot as plt
    fig, ax = plt.subplots()

    energies = np.logspace(16, 21, 101)
    ax.plot(energies, energies**3*auger_flux(energies), color='k')
    for i, (auger_comp, auger_spline, auger_label) in enumerate(
            zip(auger_components, auger_splines, auger_comp_names)):
        ax.plot(10**auger_comp[:, 0], auger_comp[:, 1], label=auger_label, color=f'C{i}', marker='o')

        ax.plot(energies, auger_spline(energies), color=f'C{i}')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('E / eV')
    ax.set_ylabel(r'$E^3 \cdot \Phi$ / (eV$^2$ km$^{-2}$ sr$^{-1}$ yr$^{-1}$)')
    ax.legend()
    fig.savefig('figs_cr/auger_composition.pdf', bbox_inches='tight')


    flux = simweights.GaisserH4a()

    fig, ax = plt.subplots()
    energies = np.logspace(5, 12, 141)
    for pdgid, _flux in zip(flux.pdgids, flux._funcs):
        ax.plot(energies, _flux(energies)*energies**3, label=pdgid.name)

    for i, (auger_comp, auger_label) in enumerate(zip(auger_components, auger_comp_names)):
        ax.plot(10**(auger_comp[:, 0]-9),
                auger_comp[:, 1] / SEC_PER_YR / 1e10 / 1e18,
                label=auger_label, color=f'C{i}', ls='--')

    ax.plot(energies,
            np.sum([_flux(energies) for _flux in flux._funcs], axis=0)*energies**3,
            color='k', label='Total')
    ax.plot(energies, (energies*1e9)**3 * auger_flux(energies*1e9) / SEC_PER_YR / 1e10 / 1e18,
            color='grey', label='Auger combined flux')
    ax.legend(bbox_to_anchor=(1,1,0,0))
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim(1e-2, 1e4)
    ax.set_xlabel('E / GeV')
    ax.set_ylabel(r'$E^3 \cdot \Phi$ / (GeV$^2$ cm$^{-2}$ sr$^{-1}$ s$^{-1}$)')
    fig.savefig('figs_cr/cosmic_ray_fluxes_GaisserH4a_auger.pdf',
                bbox_inches='tight')

    flux_ = GaisserH4aAugerZombie()
    fig, ax = plt.subplots()
    for pdgid, _flux in zip(flux_.pdgids, flux_._funcs):
        ax.plot(energies, _flux(energies)*energies**3, label=pdgid.name)

    ax.plot(energies,
            np.sum([_flux(energies) for _flux in flux._funcs], axis=0)*energies**3,
            color='k', label='Total')
    ax.plot(energies, (energies*1e9)**3 * auger_flux(energies*1e9) / SEC_PER_YR / 1e10 / 1e18,
            color='grey', label='Auger combined flux')
    ax.legend(bbox_to_anchor=(1,1,0,0))
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim(1e-2, 1e4)
    ax.set_xlabel('E / GeV')
    ax.set_ylabel(r'$E^3 \cdot \Phi$ / (GeV$^2$ cm$^{-2}$ sr$^{-1}$ s$^{-1}$)')
    fig.savefig('figs_cr/cosmic_ray_combined_fluxes_GaisserH4a_auger.pdf',
                bbox_inches='tight')


if __name__ == '__main__':
    main()
