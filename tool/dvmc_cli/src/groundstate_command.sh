# Ground state command script

run_dvmc () {
  # Constants
  N_PROC="${DVMC_MPI_PROC}"
  NAMELIST=$(find . -name "namelist.def" -type f)
  OPTIMIZED=${args[--optimized]}

  # Script
  if [ -f "${NAMELIST}" ]; then
    bold "--------------------------------------------------------------------------------"
    green_bold "[@] Static ground state numerical evaluation..."
    if [ "${OPTIMIZED}" ]; then
      if [ -f "./output/zqp_opt.dat" ]; then
        green_bold "[@] Using optimized parameters stored in './output/zqp_opt.dat'"
        mpirun -n "${N_PROC}" "${DVMC_SCRIPTS_LOCATION}"/dvmc.out "${NAMELIST}" "./output/zqp_opt.dat"
      else
        red_bold "[X] Optimized parameters: './output/zqp_opt.dat' not found in current directory."
        exit 1
      fi
    else
      mpirun -n "${N_PROC}" "${DVMC_SCRIPTS_LOCATION}"/dvmc.out "${NAMELIST}"
    fi
  else
    red_bold "[X] File named: 'namelist.def' not found in current directory."
    exit 1
  fi
}

run_dvmc
