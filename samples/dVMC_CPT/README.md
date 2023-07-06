<div align="center">

# dVMC + CPT calculation sample

</div>

This directory contains a direct application of the dVMC software usage to
compute the ground state of a small 3x4 cluster then using the Q-matrix
representation in cluster perturbation theory (CPT) to output the Fermi
surface of the system by calling PyQCM library.

## Content

- `model_3x4.py`: Defines the lattice model using PyQCM library instances.

- `fermi_surface.py`: Uses PyQCM dVMC interface to find the groundstate then
  CPT method to output the Fermi surface.

- `params`: Global dVMC solver input parameter file.

- `clean.sh`: Cleans this directory except for listed files.

- `fermi_surface_expected.pdf`: The expected Fermi surface that any user should
  obtain by running this example.

- `output/`: Working subdirectory.

## Usage

Since this example uses dVMC through PyQCM, the only thing to do is to execute
the `fermi_surface.py` script using Python 3 (with proper PyQCM installation)

#### Run script to generate A(k,w) from dVMC + CPT

```shell
python3 fermi_surface.py
```
