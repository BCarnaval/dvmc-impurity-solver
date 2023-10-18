<div align="center">

# CDMFT-dVMC calculation sample

</div>

This directory contains the minimal code to reproduce one point of 
It reproduces the results in Fig. 2 of [arxiv:2307.15738](https://arxiv.org/abs/2307.15738) using PyQCM. 
Note that in Fig. 2, the last cdmft iterations were averaged in order 
to reach better precision.

This sample is meant to be executed by a supercomputer using `run_slurm.sh` shell
script as the job submission script. The script "plot_output.py" requires the module "pandas",
which can be installed with:

```shell
pip3 pandas
```

> ### Note
>
> It takes approximately 3.5 hours to run on [Béluga](https://docs.alliancecan.ca/wiki/B%C3%A9luga/en)
> supercomputer using 32 MPI processes with 1 CPU per task and 4G of memory.

## Content

- `general_bath_1D.py`: Defines 1D lattice models using PyQCM library instances.

- `params`: Global dVMC solver input parameter file.

- `run_slurm.sh`: [SLURM](https://slurm.schedmd.com/sbatch.html) job submission bash script.

- `plot_output.py` : Python3 script to reproduce `ave_mu.pdf` from the results. Requires the "pandas" module.

- `expected/`: Directory containing expected results from the CDMFT-dVMC calculations. It also contains the ref sudirectory with the `exact_Lieb-Wu.tsv` solution.

## Usage

Since this example uses dVMC through PyQCM, the only thing to do is to execute
the `general_bath_1D.py` script using Python 3 (with proper PyQCM installation)

#### Run script to execute CDMFT-dVMC calculations

```shell
python3 general_bath_1D.py
python3 plot_output.py
```

## Post processing

To clean the directory from generated files that would be overwritten by the program,
use the `dvmc` command line interface

```shell
dvmc clean --help
```

