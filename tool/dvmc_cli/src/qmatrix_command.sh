# Q-matrix command script

q_matrix () {
  # Constants
  PROCESS="${args[processing]}"
  TOL="${args[--tolerance]}"
  USE_FILTER="${args[--use_filter]}"
  TL_FILTER="${args[--addtl_filter]}"
  K_TOL="${args[--k_tolerance]}"

  # Script selection
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

q_matrix
