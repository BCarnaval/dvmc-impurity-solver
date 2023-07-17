# Generate command script

generate () {
  # Generating global .def files
  PARAMS="${args[parameters_file]}"
  if [ -f "${PARAMS}" ]; then
    green_bold "[@] Generating .def files..."
    "${DVMC_SCRIPTS_LOCATION}"/init_params.py "${PARAMS}"
  else
    red_bold "[X] Input must be a file."
    exit 1
  fi

  # Generating excitation.def file
  if [ -z "${args[hopping_option]}" ]; then
    yellow_bold "[!] Skipping excitation.def"
  elif [ "${args[hopping_option]}" == 1 ]; then
    green_bold "[@] Using 'makeExcitation.py' to generate excitation.def"
    "${DVMC_SCRIPTS_LOCATION}"/makeExcitation.py

  elif [ "${args[hopping_option]}" == 2 ]; then
    green_bold "[@] Using 'makeExcitation_from_hopping.py' to generate excitation.def"
    "${DVMC_SCRIPTS_LOCATION}"/makeExcitation_from_hopping.py

  elif [ "${args[hopping_option]}" == 3 ]; then
    green_bold "[@] Using 'makeExcitation_from_hopping_only_t.py' to generate excitation.def"
    "${DVMC_SCRIPTS_LOCATION}"/makeExcitation_from_hopping_only_t.py
  fi
  green_bold "[@] .def files generated."
}

generate
