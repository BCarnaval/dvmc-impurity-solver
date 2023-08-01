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
    "read_gs": False, "tol": 1e-10, "save_iters": False,
    "reset_params": False, "skip_gs_iter1": False, "addtl_filter": 0,
    "pct_filter": 0.9, "cond_number": False, "use_SVD": False, "k_tol": 3,
    "const_cutoff": False, "iters_to_set_cutoff": 0, "CPT_flg": True,
    "cluster_label": 0
}


def set_dvmc_parameters(parameters: dict) -> None:
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
        version = subprocess.getoutput("dvmc --version")
        print(f"Found dVMC CLI version: {version}")
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
    verify_dvmc_installation()

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

    # Calling the dVMC solver
    iter_for_vmc = 0
    if not dvmc_parameters['exc_only'] and not dvmc_parameters['read_soln']:
        if iter_for_vmc > 0 or not dvmc_parameters['skip_gs_iter1']:
            # Ground state calculation
            if dvmc_parameters['read_gs']:
                p = subprocess.run(
                    ["dvmc", "groundstate", "--optimized"])
            else:
                p = subprocess.run(["dvmc", "groundstate"])

    if p.returncode != 0:
        print("VMC GS calculation failed")
        exit()

    # Make excitation list
    if not dvmc_parameters['gs_only'] and not dvmc_parameters['read_soln']:
        # Dynamical calculation
        p = subprocess.run(["dvmc", "excitations"])
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
        k_tol = dvmc_parameters['k_tol']
        tl_filter = dvmc_parameters['addtl_filter']
        pct_filter = dvmc_parameters['pct_filter']
        if dvmc_parameters['cond_number']:
            p = subprocess.run([
                "dvmc",
                "qmatrix",
                "cond",
                f"--tolerance={tol}",
                f"--use_filter={tl_filter}"
                f"--addtl_filter={pct_filter}"
            ])
        elif dvmc_parameters['use_SVD']:
            p = subprocess.run([
                "dvmc",
                "qmatrix",
                "svd",
                f"--tolerance={tol}",
                f"--use_filter={tl_filter}",
                f"--addtl_filter={tl_filter}",
                f"--k_tolerance={k_tol}"
            ])
        else:
            p = subprocess.run([
                "dvmc",
                "qmatrix",
                "sqrt",
                f"--tolerance={tol}",
                f"--use_filter={tl_filter}",
                f"--addtl_filter={pct_filter}"
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

    return
