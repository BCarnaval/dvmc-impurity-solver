#!/usr/bin/env bash

ITERS=$(seq 0 1 20)
PY_PATH=$(which python3)
FILE_PATH=./fermi_surface.py

main () {
  for i in ${ITERS}; do
    ntfy send QcmSimulations "Iteration ${i}..."
    ${PY_PATH} ${FILE_PATH} $i
  done

  ntfy send QcmSimulations "Simulations done!"
}

main
