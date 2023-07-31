"""Plots the Fermi surface of the model defined in 'model_3x4'
Python module using PyQCM library to invoke dVMC solver instead
of exact diagonalization.
"""
import pyqcm
import pyqcm.dvmc
import numpy as np
from model_3x4 import model
import matplotlib.pyplot as plt


def main() -> None:
    # Reading parameters as dictionnary
    params = dict(np.genfromtxt('./params', encoding='utf8', dtype=None))

    # Setting hamiltonian sector & lattice parameters
    sector = f"R0:N{params['nelec']}:S0"
    model.set_target_sectors([sector])
    model.set_parameters(
        f"""
        U={params['U']}
        t={params['t']}
        tp={params['tp']}
        tpp={params['tpp']}
        mu={params['mu']}
        """
    )

    f = open('sec', 'w')
    f.write(sector)
    f.close()

    # Setup dVMC solver through PyQCM
    pyqcm.solver = pyqcm.dvmc.dvmc_solver

    model_instance = pyqcm.model_instance(model)

    # CPT DoS
    model_instance.mdc(eta=0.12,
                       quadrant=False,
                       freq=1.5970819385867554,
                       file='fermi_surface.pdf',
                       sym='RXY'
                       )

    # Setup ED solver through PyQCM
    w, A_dvmc = model_instance.cluster_spectral_function(wmax=15, color='k')

    pyqcm.solver = None
    model_instance = pyqcm.model_instance(model)
    w, A_ed = model_instance.cluster_spectral_function(wmax=15, color='r')

    # Plotting both dVMC & ED solutions for cluster spectral functions
    dim = int(params['nelec']) // 3
    for i in range(dim):
        plt.plot(np.real(w), A_dvmc[:, i] + 2 *
                 i, color='C5', lw=2, alpha=0.95)
        plt.plot(np.real(w), A_ed[:, i] + 2 * i, color='C0')

    plt.yticks(2 * np.arange(0, dim), [str(i) for i in range(1, dim + 1)])
    plt.xlabel(r'$\omega$')
    plt.axvline(0, ls='solid', lw=0.5)
    plt.savefig("./spectrums.pdf", dpi=800)
    plt.close()

    return


if __name__ == "__main__":
    main()
