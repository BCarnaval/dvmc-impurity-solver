# Run command script

run () {
  # --------
  # Generate
  # --------
  PARAMS="${args[parameters_file]}"
  if [ -f "${PARAMS}" ]; then
    bold "--------------------------------------------------------------------------------"
    green_bold "[@] Generating .def files..."
    "${DVMC_SCRIPTS_LOCATION}"/init_params.py "${PARAMS}"
  else
    red_bold "[X] Input must be a file."
    exit 1
  fi

  # Generating excitation.def file
  green_bold "[@] Using 'makeExcitation_from_hopping_only_t.py' to generate excitation.def"
  "${DVMC_SCRIPTS_LOCATION}"/makeExcitation_from_hopping_only_t.py
  green_bold "[@] .def files generated."

  # -------------------------
  # Ground state calculations
  # -------------------------
  N_PROC="${DVMC_MPI_PROC}"
  NAMELIST=$(find . -name "namelist.def" -type f)

  # Script
  if [ -f "${NAMELIST}" ]; then
    bold "--------------------------------------------------------------------------------"
    green_bold "[@] Static ground state numerical evaluation..."
    mpirun -n "${N_PROC}" "${DVMC_SCRIPTS_LOCATION}"/dvmc.out "${NAMELIST}"
  else
    red_bold "[X] File named: 'namelist.def' not found in current directory."
    exit 1
  fi

  # ------------------------
  # Excitations calculations
  # ------------------------
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

  # --------------------------
  # Processing output binaries
  # --------------------------
  bold "--------------------------------------------------------------------------------"
  green_bold "[@] Processing ./output/zvo_nCHAm_nAHCm_0*bin files..."
  "${DVMC_SCRIPTS_LOCATION}"/mergeOutputBin.py ./output/zvo_nCHAm_nAHCm_0*bin

  # ---------------------
  # Q-matrix calculations
  # ---------------------
  PROCESS="${args[processing]}"
  TOL="${args[--tolerance]}"
  USE_FILTER="${args[--use_filter]}"
  TL_FILTER="${args[--addtl_filter]}"
  K_TOL="${args[--k_tolerance]}"

  if [ "${PROCESS}" == "cond" ]; then
    SCRIPT=dvmc_spectrum_w_cond_number.py
  elif [ "${PROCESS}" == "svd" ]; then
    SCRIPT=dvmc_spectrum_w_SVD.py
  elif [ "${PROCESS}" == "sqrt" ]; then
    SCRIPT=dvmc_spectrum_eigh_w_sqrtS.py
  else
    red_bold "[X] The option: '${PROCESS}' isn't supported. Exiting..."
    exit 1
  fi

  bold "--------------------------------------------------------------------------------"
  green_bold "[@] Computing green function using Q-matrix using script: ${SCRIPT}..."

  # Options selection
  if [[ ! "${TOL}" ]]; then
    yellow_bold "[!] '--tolerance' not specified. Using default value: 1e-10."
    TOL=1e-10
  fi
  if [[ ! "${USE_FILTER}" ]]; then
    yellow_bold "[!] '--use_filter' not specified. Using default value: 0 (false)."
    USE_FILTER=0
  fi
  if [[ ! "${TL_FILTER}" ]]; then
    yellow_bold "[!] '--addtl_filter' not specified. Using default value: 0.9."
    TL_FILTER=0.9
  fi
  if [[ "${K_TOL}" ]] && [[ ! "${PROCESS}" == "svd" ]]; then
    red_bold "[X] '--k_tolerance' flag cannot be used with process: ${PROCESS}."
    exit 1
  fi

  if [ ! "${K_TOL}" ]; then
    "${DVMC_SCRIPTS_LOCATION}"/"${SCRIPT}" spectrumpara.def output "${TOL}" "${USE_FILTER}" "${TL_FILTER}"
  else
    "${DVMC_SCRIPTS_LOCATION}"/"${SCRIPT}" spectrumpara.def output "${TOL}" "${USE_FILTER}" "${TL_FILTER}" "${K_TOL}"
  fi
  green_bold "[@] Green's function calculations done."
}

run
