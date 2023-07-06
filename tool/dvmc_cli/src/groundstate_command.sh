# Ground state command script

run_dvmc () {
  # Constants
  N_PROC="${DVMC_MPI_PROC}"
  NAMELIST="${args[namelist]}"
  OPTIMIZED="${args[optimized]}"

  # Script
  if [ $(basename "${NAMELIST}") == 'namelist.def' ] && [ -f "${NAMELIST}" ]; then
    green "[@] Static ground state numerical evaluation..."
    if [ -z "${OPTMIZED}" ]; then
      mpirun -n "${N_PROC}" "${DVMC_SCRIPTS_LOCATION}"/dvmc.out "${NAMELIST}"
    else
      green "[@] Using optimized parameters stored in '${OPTIMIZED}'"
      mpirun -n "${N_PROC}" "${DVMC_SCRIPTS_LOCATION}"/dvmc.out "${NAMELIST}" "${OPTMIZED}"
    fi
  else
    red "[X] Input must be a file named: 'namelist.def'."
    exit 1
  fi
}

run_dvmc
