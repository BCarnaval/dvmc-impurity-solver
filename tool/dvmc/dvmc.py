#!/usr/bin/env python3
import os
import pyqcm
import numpy as np


def dvmc_solver(nproc=1, gs_only=False, exc_only=False, read_soln=False,
                rand_exc=False, read_gs=True, tol=1e-5) -> None:
    """Docs
    """
    # Printing the model data in a file
    hop = pyqcm.cluster_hopping_matrix(full=True)
    np.save('hop.npy', hop)
    interactions = np.array(pyqcm.interactions())
    np.save('interaction.npy', interactions)

    a, b = pyqcm.properties()
    a = a.split()
    b = b.split()
    sec_ind = a.index('sector_1')
    sec = b[sec_ind]
    print("sector: ", sec)

    # Print parameters to write solution for QCM
    f = open('QCM_params.def', 'w')
    f.write('sector ' + sec + '\n')
    param_set = pyqcm.parameter_set()
    f.write('sector ' + sec + '\n')
    for name, val in param_set.items():
        # Print non-zero parameters
        if abs(val[0]) > 1e-8:
            if 'b' in name:
                # Print bath params without underscore
                f.write(name.partition('_')[0] + '  ' + repr(val[0]) + '\n')
            elif '_' not in name:
                # Print cluster params
                f.write(name + '  ' + repr(val[0]) + '\n')

    f.close()
    os.system("init_params.py params")

    # Calling the dvmc solver
    if not exc_only and not read_soln:
        print('Beginning ground state VMC calculation')

        # Ground state calculation
        if read_gs:
            run_cmd = "mpirun -n " + \
                repr(nproc) + " dvmc.out namelist.def ./output/zqp_opt.dat"
        else:
            run_cmd = "mpirun -n " + repr(nproc) + " dvmc.out namelist.def"

        os.system(run_cmd)

    if not gs_only and not read_soln:
        # Make excitation list
        os.system("makeExcitation_from_hopping_only_t.py")

        if rand_exc:
            os.system("randomize_excitations.py")
            os.system("mv excitation_new.def excitation.def")

        print('Beginning dynamical VMC calculation')

        # Dynamical calculation
        run_cmd = "mpirun -n " + \
            repr(nproc) + " dvmc.out namelist_G.def ./output/zqp_opt.dat"
        os.system(run_cmd)

        # Merge binary outputs
        os.system("mergeOutputBin.py output/zvo_nCHAm_nAHCm_0*bin")

        print('Calculating and printing Qmatrix')

        # Compute and output Green's functions
        os.system(
            "dvmc_spectrum_eigh_w_sqrtS.py spectrumpara.def output " + repr(tol))

    # Reading the solution
    # Here, a file qmatrix.def is read. It contains information about the
    # Hilbert space sector, the ground state energy, and the Q matrix
    if not gs_only:
        print('Reading dVMC solution')
        with open('output/qmatrix.def') as f:
            solution = f.read()
        f.close()

        pyqcm.read_cluster_model_instance(solution, 0)

    return
