# Excitations command script

excitations () {
  # Constants
  N_PROC="${DVMC_MPI_PROC}"
  OPTIMIZED="${args[optimized]}"
  NAMELIST="${args[namelist]}"

  # Script
  if [ $(basename "${NAMELIST}") == 'namelist_G.def' ] && [ -f "${NAMELIST}" ]; then
    green "[@] Beginning dynamical VMC calculations..."
    mpirun -n "${N_PROC}" "${DVMC_SCRIPTS_LOCATION}"/dvmc.out "${NAMELIST}" "${OPTIMIZED}"
  else
    red "[X] Input must be a file named: 'namelist_G.def'."
    exit 1
  fi
}

excitations
