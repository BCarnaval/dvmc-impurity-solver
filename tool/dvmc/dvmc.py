"""dVMC solver interface module for PyQCM library. Should be placed in
directory: path_to_pyqcm/qcm_wed/pyqcm/ and small changes should be made to
PyQCM __init__.py and cdmft.py modules in order to be fully functionnal.
"""
import os
import pyqcm
import subprocess
import numpy as np


# dVMC default parameters
dvmc_parameters = {
    "nproc": 1, "gs_only": False, "exc_only": False, "read_soln": False,
    "read_gs": True, "tol": 1e-10, "save_iters": False,
    "reset_params": False, "skip_gs_iter1": False, "addtl_filter": 0,
    "pct_filter": 0.9, "cond_number": False, "use_SVD": False, "k_tol": 3,
    "const_cutoff": False, "iters_to_set_cutoff": 0, "CPT_flg": True,
    "cluster_label": 0
}

# dVMC-CDMFT default parameters
dvmc_cdmft_parameters = {
    'iter_for_vmc': 0, 'E0_VMC': 0.0, 'E0_err_VMC': 0.0,
    'min_iter_E0': 5
}


def set_dvmc_parameters(parameters: dict,
                        cdmft_dvmc_params: dict = None) -> None:
    """Updates the dVMC solver parameters using input dictionnary.

    Parameters
    ----------
    parameters: dict, default=None
        User input of dVMC solver parameters.

    cdmft_params: dict, default=None
        dVMC-CDMFT specific parameters.
    """
    if pyqcm.solver != dvmc_solver:
        print("PyQCM 'solver' variable must be set to 'dvmc' "
              "in order to use this function.")
        exit()
    else:
        # Setting user's global parameters
        for key, val in parameters.items():
            if key in dvmc_parameters.keys():
                dvmc_parameters[key] = val
            else:
                pass
        # Setting user's dVMC-CDMFT parameters if given
        if cdmft_dvmc_params is not None:
            for key, val in cdmft_dvmc_params.items():
                if key in dvmc_cdmft_parameters.keys():
                    dvmc_cdmft_parameters[key] = val
                else:
                    pass
    return


def verify_dvmc_installation() -> None:
    """Determines if dVMC program is properly installed on the system.
    """
    # Calling 'dvmc' CLI to test installation
    test_help = subprocess.run(["dvmc", "--help"], stdout=subprocess.DEVNULL)
    if test_help.returncode != 0:
        print("dVMC is not properly installed on your system. Exiting...")
        exit(0)
    else:
        version = subprocess.getoutput(["dvmc --version"])
        print(f"Found dVMC CLI version: {version}")
        pyqcm.dvmc_installation = True
    return


def dvmc_solver(model_instance: pyqcm.model_instance) -> None:
    """Computes the ground state using dynamical variational Monte Carlo
    within PyQCM library, then computes the Q-matrix representation.

    Returns of the solver are internal results for PyQCM library. Specific
    dVMC results could be found inside ./output/ directory of the project.

    Parameters
    ----------
    model_instance: pyqcm.model_instance, default=None
        The PyQCM model that defines the model to solve for GS.
    """
    # Verifying dVMC installation before any further calculation
    if not pyqcm.dvmc_installation:
        verify_dvmc_installation()
    else:
        pass

    # Trying to get SLURM properties
    try:
        n_threads = os.getenv('SLURM_CPUS_PER_TASK')
        n_proc = int(os.getenv('SLURM_NTASKS'))
    except Exception:
        n_threads = 1
        n_proc = dvmc_parameters['nproc']

    # Setting ENV. variables according to threads & proc numbers
    print('nthreads = ', n_threads)
    print('nproc = ', n_proc)
    os.environ["DVMC_MPI_PROC"] = str(n_proc)
    os.environ["OMP_NUM_THREADS"] = str(n_threads)

    # Printing the model instance hopping matrix
    print('Printing hopping matrix')
    clus_lbl = dvmc_parameters['cluster_label']
    if model_instance.model.clus[clus_lbl].cluster_model.n_bath == 0:
        hop = model_instance.cluster_hopping_matrix()
    else:
        hop = model_instance.cluster_hopping_matrix(full=True)
    np.save('hop.npy', hop)

    # Printing the model instance interactions matrix
    print('Printing interaction')
    interaction = np.array(model_instance.interactions())
    np.save('interaction.npy', interaction)

    # Getting the hamiltonian sector
    print('Reading sector')
    f = open('sec', 'r')
    sec = f.readline()
    f.close()

    # Print parameters to write solution for QCM
    f = open('QCM_params.def', 'w')
    f.write('sector ' + sec + '\n')

    param_set = model_instance.model.parameter_set()
    for name, val in param_set.items():
        # Print non-zero parameters
        if abs(val[0]) > 1e-8:
            if 'b' in name:
                # Write bath params without underscore
                f.write(name.partition('_')[0]+'  '+repr(val[0])+'\n')
            elif ('_' not in name):
                # Write cluster params
                f.write(name+'  '+repr(val[0])+'\n')
    f.close()

    # .def files generation
    try:
        if not dvmc_parameters['gs_only'] and not dvmc_parameters['read_soln']:
            p = subprocess.run(["dvmc", "generate", "params", "3"])
        else:
            p = subprocess.run(["dvmc", "generate", "params"])
    except Exception:
        print("Param files not written properly")
        exit()

    if dvmc_parameters['CPT_flg']:
        iter_for_vmc = 0
    else:
        iter_for_vmc = pyqcm.cdmft.iter_for_vmc

    # Calling the dVMC solver
    if not dvmc_parameters['exc_only'] and not dvmc_parameters['read_soln']:
        if iter_for_vmc > 0 or not dvmc_parameters['skip_gs_iter1']:
            # Ground state calculation
            if dvmc_parameters['read_gs']:
                p = subprocess.run(
                    ["dvmc", "groundstate", "namelist.def", "./output/zqp_opt.dat"])
            else:
                p = subprocess.run(["dvmc", "groundstate", "namelist.def"])

    if p.returncode != 0:
        print("VMC GS calculation failed")
        exit()

    if not dvmc_parameters['CPT_flg']:
        E0_VMC, E0_err_VMC = np.loadtxt('./output/zqp_opt.dat', usecols=(0, 2))
        dvmc_cdmft_parameters['E0_VMC'] = E0_VMC
        dvmc_cdmft_parameters['E0_err_VMC'] = E0_err_VMC

    # Make excitation list
    if not dvmc_parameters['gs_only'] and not dvmc_parameters['read_soln']:
        # Dynamical calculation
        p = subprocess.run(
            ["dvmc", "excitations", "namelist_G.def", "./output/zqp_opt.dat"])
        if p.returncode != 0:
            print("Calculation of excitations failed")
            exit()

        # Merge binary outputs
        p = subprocess.run(
            ["dvmc", "process-output", "output/zvo_nCHAm_nAHCm_0"])
        if p.returncode != 0:
            print("Failed to merge binary files")
            exit()

        # Compute and output Green's functions
        tol = dvmc_parameters['tol']
        tl_filter = dvmc_parameters['addtl_filter']
        pct_filter = dvmc_parameters['pct_filter']
        if dvmc_parameters['cond_number']:
            p = subprocess.run([
                "dvmc",
                "qmatrix",
                "cond",
                repr(tol),
                repr(tl_filter),
                repr(pct_filter)
            ])
        elif dvmc_parameters['use_SVD']:
            p = subprocess.run([
                "dvmc",
                "qmatrix",
                "svd",
                repr(tol),
                repr(tl_filter),
                repr(pct_filter),
                repr(dvmc_parameters['k_tol'])
            ])
        else:
            p = subprocess.run([
                "dvmc",
                "qmatrix",
                "sqrt",
                repr(tol),
                repr(tl_filter),
                repr(pct_filter)
            ])
        if p.returncode != 0:
            print("Failed to write Q-matrix")
            exit()

    # Reading the solution here, a file qmatrix.def is read. It contains
    # information about the Hilbert space sector, the ground state energy,
    # and the Q matrix
    if not dvmc_parameters['gs_only']:
        print('Reading dVMC solution')

        with open('output/qmatrix.def') as f:
            solution = f.read()
        f.close()
        model_instance.read(solution, 0)

    print('Done reading dVMC solution')

    if dvmc_parameters['save_iters']:
        subprocess.run(["mv", "output/zqp_opt.dat",
                       "output/zqp_opt_iter"+str(iter_for_vmc)+".dat"])
        subprocess.run(["mv", "output/qmatrix.def",
                       "output/qmatrix_iter"+str(iter_for_vmc)+".def"])
        subprocess.run(["mv", "output/zbo_var_001.dat",
                       "output/zbo_var_001_iter"+str(iter_for_vmc)+".dat"])
        subprocess.run(
            ["cp", "output/S_AC.npy", "output/S_AC_"+str(iter_for_vmc)+".npy"])
        subprocess.run(
            ["cp", "output/S_CA.npy", "output/S_CA_"+str(iter_for_vmc)+".npy"])
        subprocess.run(
            ["cp", "output/H_AC.npy", "output/H_AC_"+str(iter_for_vmc)+".npy"])
        subprocess.run(
            ["cp", "output/H_CA.npy", "output/H_CA_"+str(iter_for_vmc)+".npy"])

    return
