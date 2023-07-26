# Process-output command script

process_output () {
  # Constants
  PREFIX="${args[prefix]}"

  # Script
  bold "--------------------------------------------------------------------------------"
  green_bold "[@] Processing ${PREFIX}*bin files..."
  "${DVMC_SCRIPTS_LOCATION}"/mergeOutputBin.py "${PREFIX}"*bin
}

process_output
