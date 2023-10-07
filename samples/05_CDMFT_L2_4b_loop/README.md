<div align="center">

# CDMFT-dVMC calculation sample

</div>

This directory contains the minimal code to reproduce one point of 
It reproduces the results in Fig. 2 of [arxiv:2307.15738](https://arxiv.org/abs/2307.15738) using PyQCM. 
Note that in Fig. 2, the last cdmft iterations were averaged in order 
to reach better precision.

This sample is meant to be executed by a supercomputer using `run_slurm.sh` shell
script as the job submission script. It roughly reproduces the figure 2 of the [article](https://arxiv.org/abs/2307.15738)
but with less statistics and sampling (see `expected/` directory). There is also an averaging shortcut in the
production of `./expected/ave_mu.pdf` figure compared to the original paper. Here, we used the last value outputed
from CDMFT iterations instead of avering on those iterations.

> ### Note
>
> It takes approximately 3.5 hours to run on [Béluga](https://docs.alliancecan.ca/wiki/B%C3%A9luga/en)
> supercomputer using 32 MPI processes with 1 CPU per task and 4G of memory.

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
