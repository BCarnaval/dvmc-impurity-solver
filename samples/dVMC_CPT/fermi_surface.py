"""Plots the Fermi surface of the model defined in 'model_3x4'
Python module using PyQCM library to invoke dVMC solver instead
of exact diagonalization.
"""
import pyqcm
import pyqcm.dvmc
import numpy as np
from model_3x4 import model


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

    f = open('QCM_params.def', 'w')
    f.write('sector ' + sector + '\n')
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
    return


if __name__ == "__main__":
    main()
