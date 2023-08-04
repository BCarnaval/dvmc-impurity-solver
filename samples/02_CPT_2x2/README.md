<div align="center">

# CPT-dVMC calculation sample

</div>

This directory contains a direct application of the dVMC software usage to
compute the ground state of a small 2x2 cluster then using the Q-matrix
representation in cluster perturbation theory (CPT) to output the cluster
spectral function of the system by calling PyQCM library.

This model is the same as the example contained in sample `01_dVMC`, but with
more options due to the pyqcm machinery (dos, mdc and spectrum).

## Content

- `model_2x2.py`: Defines the lattice model using PyQCM library instances.

- `cluster_spectral.py`: Uses PyQCM dVMC interface to find the groundstate then
  CPT method to output the cluster spectral function.

- `params`: Global dVMC solver input parameter file.

- `expected/`: Directory containing expected results from the CPT-dVMC calculations.

## Usage

Since this example uses dVMC through PyQCM, the only thing to do is to execute
the `cluster_spectral.py` script using Python 3 (with proper PyQCM installation)

#### Run script to generate cluster spectral function from CPT-dVMC

```shell
python3 cluster_spectral.py
```

## Post processing

To clean the directory from generated files that would be overwritten by the program,
use the `dvmc` command line interface

```shell
dvmc clean --help
```
