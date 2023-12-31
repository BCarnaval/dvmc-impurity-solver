# Clean command script

clean () {
  # Constants
  SPACE_IFS=' '
  OLD_IFS=$IFS

  if [[ "${args[spare]}" ]]; then
    yellow_bold "[!] Saved files: ${args[spare]}"

    # Reading prompt
    IFS=${SPACE_IFS}
    read -a SPARE <<< "${args[spare]}"

    # Saving those poor files...
    mkdir save
    for file in "${SPARE[@]}"; do
      if [[ -f ${file} ]]; then
        cp ${file} ./save/${file}
      else
        yellow_bold "[!] ${file} doesn't exist."
      fi
    done

    # Restoring IFS variable
    IFS=${OLD_IFS}
  fi

  # Defining trash files
  TRASH=$(find . -maxdepth 1 -type f -name "*.def" -o -name "*.dat" -o -name "*.npy" -o -name "sec" -o -name "*.pdf" -o -name "*.tsv")
  if [[ -d ./output ]]; then
    OUTPUT_TRASH=$(find ./output -type d)
  else
    OUTPUT_TRASH=''
  fi

  # Script
  bold "--------------------------------------------------------------------------------"
  if [[ -z "${TRASH}" ]]; then
    yellow_bold "[!] Nothing to clean in current directory."
  else
    green_bold "[@] Removing .def, .dat, .npy, .pdf and .tsv generated files..."
    rm -fv ${TRASH}
  fi

  if [[ -d ./save ]]; then
    for file in ./save/*; do
      NAME=$(basename ${file})
      cp ${file} ./${NAME}
    done

    rm -rf ./save
  fi

  if [[ "${args[--deep]}" ]] || [[ "${args[-d]}" ]]; then
    if [[ -z "${OUTPUT_TRASH}" ]]; then
      yellow_bold "[!] './output' directory not found."
    else
      yellow_bold "[!] Cleaning './output' directory..."
      rm -vrf ${OUTPUT_TRASH}
    fi
  fi
  green_bold "[@] Cleaning done."
}

clean
