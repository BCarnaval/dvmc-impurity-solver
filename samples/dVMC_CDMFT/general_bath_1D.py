"""Docs
"""
import pyqcm
import pyqcm.dvmc
import pyqcm.cdmft
import numpy as np


def model1D(ns: int, nb: int, S: list) -> pyqcm.lattice_model:
    """Defines a PyQCM lattice model based on given number of cluster
    sites and bath sites.

    Parameters
    ----------
    ns: int, default=None
        Number of cluster sites.

    nb: int, default=None
        Number of bath sites

    S: list, default=None
    """
    if ns % 2:
        raise ValueError('model1D must have an even number of cluster sites.')
    if nb % 2:
        raise ValueError('model1D must have an even number of bath sites.')

    # Instanciating the lattice as PyQCM cluster model
    no = ns + nb
    CM = pyqcm.cluster_model(name='clus', n_sites=ns, n_bath=nb)

    if nb != 0:
        model_name = f'{ns}L_{nb}b_gen'
        for i in range(1, nb + 1):
            CM.new_operator(op_name=f'eb{i}', op_type='one-body',
                            elem=[(i + ns, i + ns, 1),
                                  (i + ns + no, i + ns + no, 1)]
                            )
            CM.new_operator(op_name=f'tb{i}', op_type='one-body',
                            elem=[(1, i + ns, 1),
                                  (1 + no, i + ns + no, 1),
                                  (ns, i + ns, S[i - 1]),
                                  (ns + no, i + ns + no, S[i - 1])]
                            )
    else:
        model_name = f'{ns}L'

    cluster = pyqcm.cluster(X=CM, sites=[[i, 0, 0] for i in range(ns)])
    model = pyqcm.lattice_model(name=model_name,
                                clus=cluster,
                                superlattice=[[ns, 0, 0]])

    model.interaction_operator('U')
    model.hopping_operator('t', [1, 0, 0], -1)  # NN hopping

    return model


def main() -> None:
    # Reading parameters from dVMC parameter file
    params = dict(np.genfromtxt('./params', encoding='utf8', dtype=None))

    # Reading lattice properties
    ns, nb = int(params['Nc']), int(params['Nb'])
    n_electrons = ns + nb
    latt_model = model1D(ns, nb, [1, 1, -1, -1])
    latt_model.set_target_sectors([f'R0:N{n_electrons}:S0'])

    # Write sector for dVMC
    f = open('sec', 'w')
    f.write(f'R0:N{n_electrons}:S0')
    f.close()

    # Defining PyQCM parameters
    basic_params = f"""
    t = {params['t']}
    U = {params['U']}
    mu = {params['mu']}
    """

    bath_params = """
    eb1_1 = 1.3133194
    eb2_1 = -0.92430337
    eb3_1 = 1.3133423
    eb4_1 = -0.92433447
    tb1_1 = 0.42519886
    tb2_1 = 0.42648984
    tb3_1 = 0.42519751
    tb4_1 = 0.42649277
    """

    latt_model.set_parameters(basic_params + bath_params)
    varia = ['eb{:d}_1'.format(i+1) for i in range(nb)] + \
        ['tb{:d}_1'.format(i+1) for i in range(nb)]

    # dVMC parameters
    dvmc_params = {
        'save_iters': True,
        'use_SVD': True,
        'k_tol': 6
    }
    pyqcm.dvmc.dvmc_cdmft_parameters['min_iter_E0'] = 25

    # Specify PyQCM solver
    pyqcm.solver = pyqcm.dvmc.dvmc_solver
    pyqcm.dvmc.set_dvmc_parameters(parameters=dvmc_params)

    # Calling PyQCM CDMFT
    pyqcm.cdmft.CDMFT(model=latt_model, varia=varia, beta=25, accur=5e-3,
                      accur_E0=3e-3, miniter=15, maxiter=60, averages=True,
                      compute_potential_energy=True, SEF=True)
    return


if __name__ == "__main__":
    main()
