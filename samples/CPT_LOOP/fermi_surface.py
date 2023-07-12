"""Plots the Fermi surface of the model defined in 'model_3x4'
Python module using PyQCM library to invoke dVMC solver instead
of exact diagonalization.
"""
import pyqcm
import pyqcm.dvmc
import numpy as np
from model_3x4 import model


def main(iters: int) -> None:
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

    for i in range(iters):
        # Setup dVMC solver through PyQCM
        pyqcm.solver = pyqcm.dvmc.dvmc_solver
        model_instance = pyqcm.model_instance(model)

        # Cluster averages
        model_instance.properties()

        # # Cluster spectral function
        # model_instance.cluster_spectral_function(
        #     wmax=6,
        #     eta=0.12,
        #     file=f"./spectrums/spectrum_{i}.pdf")

    return


if __name__ == "__main__":
    # i = sys.argv[1]
    main(30)
