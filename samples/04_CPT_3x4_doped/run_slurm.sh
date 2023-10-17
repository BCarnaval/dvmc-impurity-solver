#!/usr/bin/env bash
#SBATCH --job-name=3x4_CPT_dVMC
#SBATCH --account=def-charleb1
#SBATCH --mem-per-cpu=1G
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --time=00-03:00           # time (DD-HH:MM)

module reset
module add python/3.11.5
module add scipy-stack/2023a
module add mkl
module add imkl
module add intel

main () {

    # Setting env variables
    export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
    export DVMC_MPI_PROC=${SLURM_NTASKS}

    # Setting working directory
    RUNDIR="${SLURM_SUBMIT_DIR}"
    echo $RUNDIR
    cd $RUNDIR

    # Execute the calculations
    python3 cluster_spectral.py
}

main
