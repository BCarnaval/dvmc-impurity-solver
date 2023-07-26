# Excitations command script

excitations () {
  # Constants
  N_PROC="${DVMC_MPI_PROC}"
  NAMELIST_G=$(find . -name "namelist_G.def" -type f)
  OPTIMIZED=$(find ./output -name "zqp_opt.dat" -type f)

  # Script
  if [[ -f "${NAMELIST_G}" ]] && [[ -f "${OPTIMIZED}" ]]; then
    bold "--------------------------------------------------------------------------------"
    green_bold "[@] Beginning dynamical VMC calculations..."
    mpirun -n "${N_PROC}" "${DVMC_SCRIPTS_LOCATION}"/dvmc.out "${NAMELIST_G}" "${OPTIMIZED}"
  elif [[ ! -f "${NAMELIST_G}" ]]; then
    red_bold "[X] File named: 'namelist_G.def' not found in current directory."
    exit 1
  elif [[ ! -f "${OPTIMIZED}" ]]; then
    red_bold "[X] File named: './output/zqp_opt.dat' not found in current directory."
    exit 1
  else
    red_bold "[X] Neither './output/zqp_opt.dat' and 'namelist_G.def' have been found."
    exit 1
  fi
}

excitations
