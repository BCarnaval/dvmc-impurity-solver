<div align="center">

# CDMFT-dVMC calculation sample

</div>

This directory contains the minimal code to reproduce one point of 
It reproduces one point of (mu,n) = (2.0,1.0) of the Fig. 2 of 
[arxiv:2307.15738](https://arxiv.org/abs/2307.15738) using PyQCM.

## Content

- `general_bath_1D.py`: Defines 1D lattice models using PyQCM library instances.

- `params`: Global dVMC solver input parameter file.

- `expected/`: Directory containing expected results from the CDMFT-dVMC calculations.

## Usage

Since this example uses dVMC through PyQCM, the only thing to do is to execute
the `general_bath_1D.py` script using Python 3 (with proper PyQCM installation)

#### Run script to execute CDMFT-dVMC calculations

```shell
python3 general_bath_1D.py
```

## Post processing

To clean the directory from generated files that would be overwritten by the program,
use the `dvmc` command line interface

```shell
dvmc clean --help
```
