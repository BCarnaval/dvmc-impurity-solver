"""Docs
"""
import os
import pyqcm
import pyqcm.dvmc
import numpy as np

from model_1D_2_4b import model


def rewrite_params(param_name: str, param_val) -> None:
    """Docs
    """
    fold = open('params', 'r')
    fnew = open('tmp_params', 'w')
    lines = fold.readlines()

    for line in lines:
        ln = line.strip('\n').split()
        if (str(param_name) in ln[0]):
            fnew.write(ln[0]+'    '+repr(param_val)+'\n')
        else:
            fnew.write(ln[0]+'    '+ln[1]+'\n')

    os.system('mv tmp_params params')
    return


np.set_printoptions(precision=4, linewidth=512, suppress=True)


def main() -> None:
    return


params = dict(np.genfromtxt('./params', encoding='utf8', dtype=None))
sector = f"R0:N{params['nelec']}:S0"
# print sector to QCM_params.def to be read by dVMC
f = open('QCM_params.def', 'w')
f.write('sector ' + sector + '\n')
f.close()

model.set_global_parameter('verbose', 0)
model.set_target_sectors([sector])
model.set_parameters(
    f"""
    U={params['U']}
    mu={params['mu']}
    t={params['t']}
    tb1_1=0.5
    tb2_1=0.5
    eb1_1=0.5
    eb2_1=-0.5
    """
)

# Instanciating the read model
model_instance = pyqcm.model_instance(model)

rewrite_params('NSROptItrStep', 3000)
rewrite_params('NSROptItrSmp', 50)

# gs_only = True
# get ground state with large
# number of optimization iterations
dvmc_solver(read_gs=False)

os.system("mv qmatrix.dat qmatrix_init.dat")
os.system("mv output/zbo_var_001.dat output/zbo_var_001_init.dat")

# set number of optimization iterations
# to smaller number, and read ground
# state from previous calculation

rewrite_params('NSROptItrStep', 1000)
rewrite_params('NSROptItrSmp', 50)
pyqcm.solver = 'dvmc'

obs_mu = observable('mu_1', 0.5e-3)
obs_t = observable('t_1', 0.5e-3)

cdmft(varia=['tb1_1', 'tb2_1', 'eb1_1', 'eb2_1'],
      accur=1e-5, accur_hybrid=1e-8, accur_dist=1e-10,
      observables=[obs_mu, obs_t])

rewrite_params('NVMCSample', 10000)
rewrite_params('NDataQtySmp', 30)
exc_only = True
# get Green's function with large
# number of samples
dvmc_solver(exc_only=exc_only)


if __name__ == "__main__":
    main()
