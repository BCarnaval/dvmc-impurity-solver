# Clean command script

clean () {
  # Constants
  TRASH=$(find . -name "*.def" -o -name "*.dat" -o -name "*.npy" -o -name "sec" -type f)
  OUTPUT_TRASH=$(find ./output -type f)

  # Script
  if [[ -z "${TRASH}" ]]; then
    yellow_bold "[!] Nothing to clean in current directory."
  else
    green_bold "[@] Removing '.def', '.dat' and '.npy' generated files..."
    rm -vf ${TRASH}
  fi

  if [[ "${args[--deep]}" ]] || [[ "${args[-d]}" ]]; then
    if [[ -z "${OUTPUT_TRASH}" ]]; then
      yellow_bold "[!] Nothing to clean in './output' directory."
      exit 0
    else
      yellow_bold "[!] Cleaning './output' directory..."
      rm -vf ${OUTPUT_TRASH}
    fi
  fi
  green_bold "[@] Cleaning done."
}

clean
