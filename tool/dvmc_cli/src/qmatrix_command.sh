# Q-matrix command script

q_matrix () {
  # Constants
  PROCESS="${args[processing]}"
  TOL="${args[tolerance]}"
  USE_FILTER="${args[use_filter]}"
  TL_FILTER="${args[addtl_filter]}"
  K_TOL="${args[k_tolerance]}"

  # Script
  if [ "${PROCESS}" == "cond" ]; then
    SCRIPT=dvmc_spectrum_w_cond_number.py
  elif [ "${PROCESS}" == "svd" ]; then
    SCRIPT=dvmc_spectrum_w_SVD.py
  elif [ "${PROCESS}" == "sqrt" ]; then
    SCRIPT=dvmc_spectrum_eigh_w_sqrtS.py
  elif [ -z "${PROCESS}" ]; then
    SCRIPT=dvmc_spectrum.py
  else
    red "[X] The option: '${PROCESS}' isn't supported. Exiting..."
    exit 1
  fi

  green "[@] Computing green function using Q-matrix using script: ${SCRIPT}..."
  if [ -z "${K_TOL}" ]; then
    "${DVMC_SCRIPTS_LOCATION}"/"${SCRIPT}" spectrumpara.def output "${TOL}" "${USE_FILTER}" "${TL_FILTER}"
  else
    "${DVMC_SCRIPTS_LOCATION}"/"${SCRIPT}" spectrumpara.def output "${TOL}" "${USE_FILTER}" "${TL_FILTER}" "${K_TOL}"
  fi
  green "[@] Green's function calculations done."
}

q_matrix
