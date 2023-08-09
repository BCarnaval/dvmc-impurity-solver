"""Plots the cluster spectral function of the model defined in 'model_3x4'
Python module using PyQCM library to invoke dVMC solver instead
of exact diagonalization.
"""
import pyqcm
import numpy as np
from model_3x4 import model
from dvmc import dvmc_solver
import matplotlib.pyplot as plt


def compute_spectral_weight(model_inst: pyqcm.model_instance,
                            wmax: float = 10.0,
                            eta: float = 0.05) -> tuple[np.ndarray]:
    """Locally computing spectral weight for given PyQCM
    model instance.

    Parameters
    ----------
    model_inst: pyqcm.model_instance, default=None
        Model for which compute the spectral weight.

    wmax: float, default=10.0
        Range of frequencies treated as (-wmax, wmax).

    eta: float, default=0.1
        Lorentzian broadening.

    Returns
    -------
    A: np.ndarray
        Spectral weight.

    w: np.ndarray
        Frequencies array.
    """
    clus = 0
    d = model_inst.model.dimGFC[clus]

    # Computing frequencies real & complex arrays
    w = np.arange(-wmax, wmax + 1e-6, eta / 4.0)
    wc = np.array([x + eta*1j for x in w], dtype=complex)

    # Computing spectral weight using cluster Green function
    A = np.zeros((len(w), d))
    for i in range(len(w)):
        g = model_inst.cluster_Green_function(wc[i], clus, False, False)
        for j in range(d):
            A[i, j] += -g[j, j].imag

    return w, A


def plot_results() -> None:
    """Plots the cluster spectral functions for both ED and
    dVMC impurity solver.
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
    w_dvmc, A_dvmc = compute_spectral_weight(model_dvmc, 15.0, 0.1)

    # ED part
    model_ed = pyqcm.model_instance(model)
    w_ed, A_ed = compute_spectral_weight(model_ed, 15.0)

    # Plotting both dVMC & ED solutions for cluster spectral functions
    dim = int(params['nelec']) // 3
    for i in range(dim):
        plt.plot(np.real(w_dvmc), A_dvmc[:, i] + 1.5 *
                 i, lw=2.5, alpha=1,
                 label='dVMC' if i == 0 else '_nolabel_')
        plt.plot(np.real(w_ed), A_ed[:, i] + 1.5 * i,
                 color='k', label='ED' if i == 0 else '_nolabel_')

    # V/H lines setup
    [plt.axhline(y=i, ls='solid', color='k')
     for i in 1.5 * np.array([0, 1, 2, 3])]
    plt.axvline(0, ls='dashed', lw=0.5, color='k')

    plt.yticks(1.5 * np.arange(0, dim), [str(i) for i in range(1, dim + 1)])
    plt.xlabel(r'$\omega$')
    plt.legend()
    plt.savefig("./spectrums.pdf", dpi=800)
    plt.close()

    return


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
    pyqcm.model_instance(model)

    plot_results()

    return


if __name__ == "__main__":
    main()
