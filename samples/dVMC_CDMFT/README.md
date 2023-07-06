<div align="center">

# dVMC + CDMFT calculation sample

</div>

This directory contains a direct application of the dVMC software usage to
compute the ground state of a small 1D cluster with bath sites then using
the Q-matrix representation in cluster dynamical mean field theory (CDMFT)
from PyQCM library.

## Content

- `general_bath_1D.py`: Defines 1D lattice models using PyQCM library instances.

- `params`: Global dVMC solver input parameter file.

- `clean.sh`: Cleans this directory except for listed files.

- `spectrum_rspace_expected.pdf`: The expected spectrum that any user should
  obtain by running this example.

- `output/`: Working subdirectory.

## Usage

Since this example uses dVMC through PyQCM, the only thing to do is to execute
the `general_bath_1D.py` script using Python 3 (with proper PyQCM installation)

#### Run script to generate A(k,w) from dVMC + CDMFT

```shell
python3 general_bath_1D.py
```
