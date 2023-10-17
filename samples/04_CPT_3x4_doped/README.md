
<div align="center">

# CPT-dVMC calculation sample

</div>

This directory contains the minimal code to reproduce the results of
Fig. 3(a) of [arxiv:2209.08092v1](https://arxiv.org/abs/2209.08092v1) (Fig. 8(a) of [PhysRevB.106.245132](https://doi.org/10.1103/PhysRevB.106.245132)) using PyQCM.

It computes the Fermi surface of a small hole doped 3x4 cluster using 
perturbation theory (CPT).

> ### Note
>
> It takes approximately 40 minutes to run [Béluga](https://docs.alliancecan.ca/wiki/B%C3%A9luga/en)
> supercomputer using 64 MPI processes with 1 CPU per task and 1G of memory.

## Content

- `model_3x4.py`: Defines the lattice model using PyQCM library instances.

- `cluster_spectral.py`: Uses PyQCM dVMC interface to find the groundstate then
  CPT method to output the cluster spectral function of the system.

- `params`: Global dVMC solver input parameter file.

- `expected/`: Directory containing expected results from the CPT-dVMC calculations.

- `run_slurm.sh`: [SLURM](https://slurm.schedmd.com/sbatch.html) job submission bash script.

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
