<div align="center">

# CDMFT-dVMC calculation sample

</div>

This directory contains a direct application of the dVMC software to
compute the ground state of a 1D cluster of 2 sites with 4 bath sites then using
the Q-matrix representation in cluster dynamical mean field theory (CDMFT)
from PyQCM library.

This sample is meant to be executed by a supercomputer using `run_slurm.sh` shell
script as the job submission script. It roughly reproduces the figure 2 of the article
but with less statistics and sampling (see `expected/` directory).

> ### Note
>
> It takes approximately 3.5 hours to run on [Béluga](https://docs.alliancecan.ca/wiki/B%C3%A9luga/en)
> supercomputer.

## Content

- `general_bath_1D.py`: Defines 1D lattice models using PyQCM library instances.

- `params`: Global dVMC solver input parameter file.

- `run_slurm.sh`: [SLURM](https://slurm.schedmd.com/sbatch.html) job submission bash script.

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
