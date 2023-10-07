"""Plots the cluster spectral function of the model defined in 'model_2x2'
Python module using PyQCM library to invoke dVMC solver instead
of exact diagonalization.
"""
import pyqcm
import numpy as np
from model_2x2 import model
from dvmc import dvmc_solver


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
    pyqcm.solver = dvmc_solver
    model_instance = pyqcm.model_instance(model)

    # Computing cluster spectral functions
    model_instance.cluster_spectral_function(
        wmax=15, eta=0.1, file='cluster_spectral.pdf')

    # Computing Fermi surface
    model_instance.mdc(nk=200, file='fermi_surface.pdf', sym='RXY')

    return


if __name__ == "__main__":
    main()
