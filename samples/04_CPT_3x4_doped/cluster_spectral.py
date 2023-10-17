"""Plots the Fermi surface of the model defined in 'model_3x4'
Python module using PyQCM library to invoke dVMC solver instead
of exact diagonalization.
"""
import pyqcm
import numpy as np
from model_3x4 import model
from dvmc import dvmc_solver
import matplotlib.pyplot as plt
from pyqcm import *

def plot_FS() -> None:
    """Computes the Fermi energy and plots the Fermi surface.
    """
    # Reading parameters as dictionnary
    pyqcm.solver = None
    params = dict(np.genfromtxt('./params', encoding='utf8', dtype=None))

    # dVMC part
    model_dvmc = pyqcm.model_instance(model)
    with open('./output/qmatrix.def') as qmatrix:
        dvmc_sol = qmatrix.read()
        qmatrix.close()

    model_dvmc.read(dvmc_sol, 0)
    #w_dvmc, A_dvmc = compute_spectral_weight(model_dvmc,15.0, 0.1)

    # compute Fermi energy
    Ne = int(params['nelec'])
    N  = int(params['L'])*int(params['W'])

    w = np.linspace(-30,30,1000)
    model_dvmc.plot_DoS(w=w)

    data = np.loadtxt('dos.tsv',skiprows=1)
    frequencies    = data[:,0]
    cumulative_dos = data[:,2]
    
    # find the chemical potential corresponding to this specific mu.
    mu = np.interp(0.5*(float(Ne)/float(N)),cumulative_dos,frequencies)

    # Computing Fermi surface
    model_dvmc.mdc(nk=100, file='fermi_surface.pdf', sym='RXY',freq=mu)
        
    return


def main() -> None:
    
    skip_calculations = True
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
    if not skip_calculations:  # we assume that the solution has been obtained before
        pyqcm.solver = dvmc_solver
        model_instance = pyqcm.model_instance(model)
    
    plot_FS()

    return


if __name__ == "__main__":
    main()
