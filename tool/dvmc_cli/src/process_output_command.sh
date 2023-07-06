# Process-output command script

process_output () {
  PREFIX="${args[prefix]}"
  green "[@] Processing ${PREFIX}*bin files..."
  "${DVMC_SCRIPTS_LOCATION}"/mergeOutputBin.py "${PREFIX}"*bin
}

process_output
