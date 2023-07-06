#!/usr/bin/env bash
#
# Example output files cleaning bash script
#
# Author: Antoine de Lagrave (BCarnaval)

RED="$(printf '\033[31m')"
CYAN="$(printf '\033[36m')"
GREEN="$(printf '\033[32m')"
WHITE="$(printf '\033[37m')"
ORANGE="$(printf '\033[33m')"

# Turn on extended globbing in bash shell
shopt -s extglob

clean_inputs () {
  echo -e "${WHITE}-----------------------------------------------"
  echo -e "${GREEN}[@] Deleting '.dat', '.pdf' and '.def' files..."
  rm -v !(model_3x4.py|README.md|params|fermi_surface.py|clean.sh|*.pdf)
  echo -e "${GREEN}[@] Input files deleted."
  echo -e "${WHITE}-----------------------------------------------"
}

clean_output () {
  echo -e "${GREEN}[@] Deleting './output/ directory content...'"
  rm ./output/*
  rm -rf /__pycache__/
  echo -e "${GREEN}[@] './output/' directory emptied."
  echo -e "${WHITE}-----------------------------------------------"
}

main () {
  clean_inputs
  clean_output
}

main
