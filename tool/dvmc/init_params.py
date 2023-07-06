#!/usr/bin/env python3

import os
import sys
import hopping
import numpy as np
import numpy.linalg as la

"""
opt_fij = {'all', 'cluster_bath', 'cluster_only'}
opt_jastrow = {'all', 'cluster_bath', 'cluster_only'}
opt_gutz = {'all', 'cluster_only'}

'all' --> optimize all parameters
'cluster_bath' --> optimize all cluster-cluster and cluster-bath parameters
'cluster_only' --> optimize only cluster-cluster parameters

input_fij = {0, 1}

0 --> use default initialization for f_ij
1 --> compute f_ij by diagonalizing non-interacting Hamiltionian
(see Eq.5.19 from mVMC user guide)
"""

param_list = [
        'Ndim', 'L', 'W', 'Nb', 'Nc', 'Lsub', 'Wsub', 'U',
        't', 'tp', 'tpp', 'mu', 'tb', 'eb', 'nelec', 'readhop',
        'input_fij', 'random_fij', 'break_sym', 'pinning_AFM',
        'treat_fij_open_shell', 'fij_NI', 'delta_dens', 'BCs',
        'read_fij', 'CDataFileHead', 'CParaFileHead', 'NVMCCalMode',
        'NLanczosMode', 'NDataIdxStart', 'NDataQtySmp',
        'NSPGaussLeg', 'NMPTrans', 'NSROptItrStep', 'NSROptItrSmp',
        'DSROptRedCut', 'DSROptStaDel', 'DSROptStepDt', 'NVMCWarmUp',
        'NVMCInterval', 'NVMCSample', 'NExUpdatePath', 'RndSeed',
        'NSplitSize', 'NStore', 'NSRCG', '2Sz', 'NSPStot',
        'add_noise_fij', 'run_dvmc', 'eta', 'w_min', 'w_max',
        'w_min_data', 'w_max_data', 'Nw', 'exc_choice', 'kPath', 'nh',
        'exc_per_site', 'sym_fij', 'opt_fij', 'opt_jastrow',
        'opt_gutz', 'use_jastrow', 'avg_prev_iters', 'alpha', 'beta'
        ]

default_vals = [
        2, 1, 1, 0, 1, 1, 1, 0.0,
        1.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 2, 1, 0, 0, 0, 0.0,
        0, 1, 0.2, 'periodic',
        0, 'zbo', 'zqp', 0, 0, 1,
        20, 8, 1, 100, 10, 1e-5, 0.01,
        0.01, 100, 1, 200, 0,
        123456789, 1, 1, 0, 0, 0, 0, 0,
        0.1, -8.0, 12.0, -15.0, 15.0, 2000,
        'all', '0', 1, 100, 0, 'all',
        'all', 'all', True, 0, 0.5, 0.5
        ]


def read_params(input_file):
    """Docs
    """
    not_read = np.zeros(len(param_list), dtype=int)
    not_read.fill(1)
    out_params = int(len(param_list))*[0]
    f = open(input_file, 'r')
    for line in f:
        param_name = line.split()[0]
        param_val = line.split()[1]
        for param in param_list:
            if param_name == param:
                ind = param_list.index(param)
                out_params[ind] = param_val
                not_read[ind] = 0
    for param in param_list:
        ind = param_list.index(param)
        if not_read[ind]:
            out_params[ind] = default_vals[ind]
    return out_params


params = read_params(sys.argv[1])

# Cluster dimensions
Ndim = int(params[param_list.index('Ndim')])
Lx = int(params[param_list.index('W')])
Ly = int(params[param_list.index('L')])

# Number of bath sites
Nb = int(params[param_list.index('Nb')])
Nc = int(params[param_list.index('Nc')])
Ns = Nc + Nb
print("Nb=%d Nc=%d Ns=%d" % (Nb, Nc, Ns))

# Interaction strength
Uhubb = float(params[param_list.index('U')])

readhop = int(params[param_list.index('readhop')])
BCs = params[param_list.index('BCs')]

# If hopping is not uniform, read hopping elements
avg_prev_iters = int(params[param_list.index('avg_prev_iters')])

if avg_prev_iters:
    alpha = float(params[param_list.index('alpha')])
    beta = float(params[param_list.index('beta')])

v0 = float(params[param_list.index('pinning_AFM')])

if readhop == 1:
    H0 = np.load('hop.npy')
    if avg_prev_iters and os.path.exists('hop_prev_iter.npy'):
        print("averaging previous hopping matrices")
        print(H0)
        H0_prev = np.load('hop_prev_iter.npy')
        print(H0_prev)
        H0 = alpha*H0+beta*H0_prev
        print(H0)

    # Difference in sign convention between qcm and mVMC
    H0 = -1.0*H0

# Hopping, chem potential, bath energy
if readhop != 1:
    t = float(params[param_list.index('t')])
    tp = float(params[param_list.index('tp')])
    tpp = float(params[param_list.index('tpp')])
    mu = float(params[param_list.index('mu')])
    print('t   = ', t)
    print('tp  = ', tp)
    print('tpp = ', tpp)

    if Nb > 0:
        bath_params = np.loadtxt('bath_params')
        tb_inds = np.where(bath_params[:, 0] != bath_params[:, 1])[0]
        eb_inds = np.where(bath_params[:, 0] == bath_params[:, 1])[0]
        tb = bath_params[tb_inds]
        eb = bath_params[eb_inds]
    else:
        tb = [0, 0]
        eb = [0, 0]

    # Set hopping
    if 'per' in BCs:
        H0 = hopping.set_hopping(t, mu, tb, eb, Nc, Nb, Lx, Ly)
    elif 'ope' in BCs:
        H0 = hopping.set_hopping_open(t, mu, tb, eb, Nc, Nb, Lx, Ly, tp=tp, tpp=tpp)
    else:
        print("Boundary condtions not set properly!")
        sys.exit()

    v, w = la.eigh(H0, UPLO='U')
    for ev in v:
        print(ev)

    print(2*np.sum(v[:4]))
    Vpin = np.zeros(H0.shape[0])
    if abs(v0) > 1e-10:
        print('applying pinning field of strength :', v0)
        Vpin[:Lx:2] = -1
        Vpin[1:Lx:2] = 1
        Vpin *= v0

# Write hop.dat
fh = open('hop.dat', 'w')
for sp in range(2):
    for i in range(H0.shape[0]):
        for j in range(H0.shape[1]):
            fh.write('%d  %d  %d  %d  % 6.4f  % 6.4f\n' % (
                i, sp, j, sp, H0[i, j].real, H0[i, j].imag))

fh.close()

break_sym = int(params[param_list.index('break_sym')])

# Add small symmetry-breaking term (useful for open shells)
if break_sym:
    for i in range(H0.shape[0]):
        H0[i, i] += 1e-6 * (-1)**i

# Number of electrons
ne = int(params[param_list.index('nelec')])

# Whether or not to set initial fij manually
input_fij = int(params[param_list.index('input_fij')])

# Whether or not to read initial fij
read_fij = int(params[param_list.index('read_fij')])

if input_fij and read_fij:
    print("Input fij and read fij cannot both be true")

# Whether or not to symmetrize f_ij
sym_fij = int(params[param_list.index('sym_fij')])
fij_NI = int(params[param_list.index('fij_NI')])

Lsub = int(params[param_list.index('Lsub')])
Wsub = int(params[param_list.index('Wsub')])
Nsub = Lsub * Wsub
ntrans = Lsub * Wsub

# Keep list of site numbers in 1st sublattice
sublatt_sitenums = []

# Keep list of site numbers outside 1st sublattice
remain_sitenums = []

# Coords of lattice sites
coords = np.zeros((Nc, Ndim), dtype=int)
for i in range(Nc):
    # hop.coor(Nc, Ndim, i, [Lx, Ly])
    coords[i, :] = hopping.coor(Nc, Ndim, i, [Lx, Ly])
    if coords[i, 0] < Wsub and coords[i, 1] < Lsub:
        sublatt_sitenums.append(i)
    else:
        remain_sitenums.append(i)

# Append bath sites to list of site numbers outside 1st sublattice
bath_sitenums = range(Nc, Ns)
remain_sitenums.extend(bath_sitenums)

# Write files for first VMC calculation
orbitalidx = open('orbitalidx.def', 'w')

# Symmetrize cluster-bath
cluster_bath_sym = False

# Symmetrize bath-bath
bath_bath_sym = False

# Whether to optimize bath-bath (cluster-bath is always optimized)
bath_bath_opt = False

# Minimum number of orbitals
norbs = Ns**2

fij_ind = np.zeros((Ns, Ns), dtype=int)

# Fill f_ij for i in first sublattice
ind = 0
for i in sublatt_sitenums:
    for j in range(Ns):
        fij_ind[i, j] = ind
        ind += 1

for i in remain_sitenums:
    for j in range(Ns):
        if i < Nc and j < Nc:
            # Symmetrize cluster-cluster
            fold_x, fold_y = hopping.fold_vector(coords[i, 0], coords[i, 1], Wsub, Lsub)
            sitej_folded = hopping.fold(coords[j, 0], coords[j, 1], fold_x, fold_y, Lx, Ly)
            sitei_folded = hopping.fold(coords[i, 0], coords[i, 1], fold_x, fold_y, Lx, Ly)
            tmp = fij_ind[sitei_folded, sitej_folded]
            fij_ind[i, j] = tmp

        elif i < Nc and j >= Nc:
            if cluster_bath_sym:
                # Symmetrize cluster-bath
                sitei_folded = hopping.fold(coords[i, 0], coords[i, 1], fold_x, fold_y, Lx, Ly)
                tmp = fij_ind[sitei_folded, j]
                fij_ind[i, j] = tmp
            else:
                fij_ind[i, j] = ind
                ind += 1

        if i >= Nc and j < Nc:
            if cluster_bath_sym:
                # Symmetrize cluster-bath
                sitej_folded = hopping.fold(coords[j, 0], coords[j, 1], fold_x, fold_y, Lx, Ly)
                tmp = fij_ind[i, sitej_folded]
                fij_ind[i, j] = tmp
            else:
                fij_ind[i, j] = ind
                ind += 1

        elif i >= Nc and j >= Nc:
            if bath_bath_sym:
                print("Bath-bath symmetry of fij not supported")
                sys.exit()
            else:
                fij_ind[i, j] = ind
                ind += 1

if sym_fij:
    import copy
    import fij_sym
    fij_copy = copy.deepcopy(fij_ind)
    fij_ind, norbs = fij_sym.sym()

# Complex type (1 --> complex, 0 --> real)
ctype = 0
orbitalidx.write('=============================================\n')
orbitalidx.write('NorbitalIdx\t\t'+repr(norbs)+'\n')
orbitalidx.write('ComplexType\t\t'+repr(ctype)+'\n')
orbitalidx.write('=============================================\n')
orbitalidx.write('=============================================\n')

# Which parameters to optimize (1 --> optimize, 0 --> do not optimize)
opt_flg = np.zeros(norbs, dtype=int)
opt_fij = repr(params[param_list.index('opt_fij')])
ind = 0
for i in range(Ns):
    for j in range(Ns):
        orbitalidx.write('\t'+repr(i)+'\t'+repr(j)+'\t'+repr(fij_ind[i, j])+'\n')
        if 'all' in opt_fij:
            opt_flg[fij_ind[i, j]] = 1

        elif 'cluster_bath' in opt_fij:
            if i >= Nc and j >= Nc:
                opt_flg[fij_ind[i, j]] = 0
            else:
                opt_flg[fij_ind[i, j]] = 1

        elif 'cluster_only' in opt_fij:
            if i >= Nc or j >= Nc:
                opt_flg[fij_ind[i, j]] = 0
            else:
                opt_flg[fij_ind[i, j]] = 1

        else:
            print("opt_fij flag must be 'all', 'cluster_bath', or 'cluster_only'")
            sys.exit()

for i in range(norbs):
    orbitalidx.write('\t'+repr(i)+'\t'+repr(opt_flg[i])+'\n')
orbitalidx.close()

locspn = open('locspn.def', 'w')

# Number of pinned spins (set to 0 for now)
nlocspin = 0
if nlocspin != 0:
    print('Localized spins not allowed at the moment')
    sys.exit()
else:
    locspin_site = np.zeros(Ns, dtype=int)
    locspin_site.fill(0)

# Write locspn.def
locspn.write('================================\n')
locspn.write('Nlocspn\t\t'+repr(nlocspin)+'\n')
locspn.write('================================\n')
locspn.write('========i_0LocSpn_1IteElc ======\n')
locspn.write('================================\n')
for i in range(Ns):
    locspn.write('\t'+repr(i)+'\t'+repr(locspin_site[i])+'\n')
locspn.close()


coulombintra = open('coulombintra.def', 'w')

# Write coulombintra.def
coulombintra.write('=============================================\n')
coulombintra.write('NCoulombIntra\t\t'+repr(Ns)+'\n')
coulombintra.write('=============================================\n')
coulombintra.write('================ CoulombIntra ===============\n')
coulombintra.write('=============================================\n')

if readhop:
    Uc = np.load('interaction.npy')
else:
    Uc = np.zeros(Ns)
    Uc[:Nc].fill(Uhubb)

print(Uc)

if readhop:
    for i in range(Nc):
        coulombintra.write('\t' + repr(int(Uc[i, 0])) + '\t' + repr(Uc[i, 2]) + '\n')
else:
    for i in range(Nc):
        coulombintra.write('\t' + repr(int(i)) + '\t' + repr(Uc[i]) + '\n')
for i in range(Nc, Ns, 1):
    coulombintra.write('\t' + repr(i) + '\t' + '0' + '\n')
coulombintra.close()

qptransidx = open('qptransidx.def', 'w')
ntrans = 1

# Weights of symmetry transformations
w_trans = np.zeros(ntrans)

# Set all weights to 1.0 for now
w_trans.fill(1.0)

# Write qptransidx.def
qptransidx.write('=============================================\n')
qptransidx.write('NQPTrans' + '\t\t' + repr(ntrans) + '\n')
qptransidx.write('=============================================\n')
qptransidx.write('======== TrIdx_TrWeight_and_TrIdx_i_xi ======\n')
qptransidx.write('=============================================\n')

for i in range(ntrans):
    qptransidx.write('\t' + repr(i) + '\t' + repr(1.0) + '\n')

for i in range(Ns):
    qptransidx.write('\t%d\t%d\t%d\n' % (0, i, i))

qptransidx.close()

trans = open('trans.def', 'w')
ntransfer = 0  # 2*Ns*Ns

for i in range(H0.shape[0]):
    for j in range(H0.shape[1]):
        if i == j:
            ntransfer += 1
        elif np.abs(H0[i, j]) > 1e-6:
            ntransfer += 1

# Factor of 2 for spin up and spin down
ntransfer *= 2

# Write trans.def
trans.write('========================\n')
trans.write('NTransfer' + '\t' + repr(ntransfer) + '\n')
trans.write('========================\n')
trans.write('========i_j_s_tijs======\n')
trans.write('========================\n')

H0_orig = H0
# Loop over spins (currently only spin-conserved hopping)
for s in range(2):
    if np.abs(v0) > 1e-9:
        H0 = H0_orig + np.diag((-1)**s*Vpin)

    for i in range(Ns):
        for j in range(Ns):
            if i == j:
                trans.write('\t' + repr(i) + '\t' + repr(s) + '\t' + repr(j) + '\t' + repr(s) + '\t' + repr(H0[i, j].real) + '\t'
                            + repr(H0[i, j].imag) + '\n')
            elif np.abs(H0[i, j]) > 1e-6:
                trans.write('\t' + repr(i) + '\t' + repr(s) + '\t' + repr(j) + '\t' + repr(s) + '\t' + repr(H0[i, j].real) + '\t'
                            + repr(H0[i, j].imag) + '\n')
trans.close()

gutz = open('gutzwilleridx.def', 'w')

# Number of gutzwiller params = number of sublattice sites + number of
# bath sites
ngutz = Ns

if sym_fij:
    import gi_sym
    gi_ind, ngutz = gi_sym.sym()

# Write gutzwiller.def
gutz.write('===========================\n')
gutz.write('NGutzwillerIdx' + '\t' + repr(ngutz) + '\n')
gutz.write('ComplexType' + '\t' + repr(0) + '\n')
gutz.write('===========================\n')
gutz.write('===========================\n')

for i in range(Ns):
    gutz.write('\t%d\t%d\n' % (i, i))
for i in range(Ns):
    gutz.write('\t%d\t%d\n' % (i, 1))

# Optimize gutzwiller params on cluster
opt = [0]*ngutz
opt_gutz = repr(params[param_list.index('opt_gutz')])
gutz.close()

jastrow = open('jastrowidx.def', 'w')

# Number of jastrow params
jastrowidx = np.zeros((Ns, Ns), dtype=int)

# Copy orbital index
for i in range(Ns):
    for j in range(Ns):
        jastrowidx[i, j] = fij_ind[i, j]

# Symmetrize v_ij = v_ji
for iorb in range(norbs):
    for i in range(Ns):
        for j in range(Ns):
            if jastrowidx[i, j] == iorb:
                jastrowidx[j, i] = jastrowidx[i, j]

njastrow = 0
for isite in range(Ns):
    for jsite in range(isite):
        if (jastrowidx[isite, jsite] >= 0):
            iJastrow = jastrowidx[isite, jsite]
            njastrow -= 1
            for isite1 in range(Ns):
                for jsite1 in range(Ns):
                    if (jastrowidx[isite1, jsite1] == iJastrow):
                        jastrowidx[isite1, jsite1] = njastrow


njastrow = -njastrow
opt_flg = np.zeros(njastrow, dtype=int)
opt_jastrow = repr(params[param_list.index('opt_jastrow')])
use_jastrow = repr(params[param_list.index('use_jastrow')])
if 'false' in use_jastrow or 'False' in use_jastrow:
    use_jastrow = False

ind = 0
for isite in range(Ns):
    for jsite in range(Ns):
        jastrowidx[isite, jsite] = -1 - jastrowidx[isite, jsite]
        if isite != jsite:
            if 'all' in opt_jastrow:
                opt_flg[jastrowidx[isite, jsite]] = 1

            elif 'cluster_bath' in opt_jastrow:
                if isite >= Nc and jsite >= Nc:
                    opt_flg[jastrowidx[isite, jsite]] = 0
                else:
                    opt_flg[jastrowidx[isite, jsite]] = 1

            elif 'cluster_only' in opt_jastrow:
                if isite >= Nc or jsite >= Nc:
                    opt_flg[jastrowidx[isite, jsite]] = 0
                else:
                    opt_flg[jastrowidx[isite, jsite]] = 1

            else:
                print("opt_jastrow flag must be 'all', 'cluster_bath', or 'cluster_only'")
                sys.exit()


# Append bath jastrow labels (these will not be optimized)
njastrow = njastrow

# Write jastrow.def
jastrow.write('===========================\n')
jastrow.write('NJastrowIdx'+'\t'+repr(njastrow)+'\n')
jastrow.write('ComplexType'+'\t'+repr(0)+'\n')
jastrow.write('===========================\n')
jastrow.write('===========================\n')

for i in range(Ns):
    for j in range(Ns):
        if i != j:
            jastrow.write('\t'+repr(i)+'\t'+repr(j)+'\t'+repr(jastrowidx[i, j])+'\n')

for i in range(njastrow):
    jastrow.write('\t'+repr(i)+'\t'+repr(opt_flg[i])+'\n')
    ind += 1

jastrow.close()

if input_fij:
    orbital = open('orbital.def', 'w')

    random_fij = int(params[param_list.index('random_fij')])
    if random_fij:
        scale = 0.2
        from numpy.random import default_rng
        rng = default_rng()
        fij = np.zeros(norbs, dtype=complex)

        # Random floats on interal scale*[-1,1)
        a = scale * 2 * (rng.random(norbs) - 0.5)
        fij = a
    else:
        treat_open_shell = int(params[param_list.index('treat_fij_open_shell')])

        add_noise_fij = int(params[param_list.index('add_noise_fij')])
        if fij_NI:
            fij = hopping.fij_from_H0(-1.0*H0, ne, treat_open_shell)
        else:
            avg_dens = ne/Ns
            delta_dens = float(params[param_list.index('delta_dens')])
            fij = hopping.fij_from_H_MF(-1.0*H0, Lx, Ly, ne, Uhubb, avg_dens, delta_dens)

        if add_noise_fij:
            dn = 0.1
            fij = fij + dn*np.random.random(len(fij))

    # Write orbital.def
    orbital.write('========================\n')
    orbital.write('NOrbitalIDx' + '\t' + repr(norbs) + '\n')
    orbital.write('========================\n')
    orbital.write('========i_j_s_tijs======\n')
    orbital.write('========================\n')

    fij_ind_flat = fij_ind.flatten()
    if sym_fij:
        orb_list = []
        fij_new = np.zeros(norbs, dtype=complex)
        ii = 0
        for i in range(Ns**2):
            print(i, fij[i].real)
            if(fij_ind_flat[i] not in orb_list):
                fij_new[fij_ind_flat[i]] = fij[i]
                orb_list.append(fij_ind_flat[i])
        fij = fij_new
    for i in range(norbs):
        orbital.write(repr(i)+'\t'+repr(fij[i].real)+'\t'+repr(fij[i].imag)+'\n')
    orbital.close()

modpara = open('modpara.def', 'w')

datafile_head = str(params[param_list.index('CDataFileHead')])
parafile_head = str(params[param_list.index('CParaFileHead')])
calc_mode = int(params[param_list.index('NVMCCalMode')])
lanczos_mode = int(params[param_list.index('NLanczosMode')])
nidx_start = int(params[param_list.index('NDataIdxStart')])
nsamp = int(params[param_list.index('NDataQtySmp')])
Sz = int(params[param_list.index('2Sz')])
ngauss = int(params[param_list.index('NSPGaussLeg')])
nspstot = int(params[param_list.index('NSPStot')])
nmptrans = int(params[param_list.index('NMPTrans')])
nSRsteps = int(params[param_list.index('NSROptItrStep')])
nSRsamp = int(params[param_list.index('NSROptItrSmp')])
SRredcut = float(params[param_list.index('DSROptRedCut')])
SRstadel = float(params[param_list.index('DSROptStaDel')])
SRdt = float(params[param_list.index('DSROptStepDt')])
nwarmup = int(params[param_list.index('NVMCWarmUp')])
ninterval = int(params[param_list.index('NVMCInterval')])
nVMCsamp = int(params[param_list.index('NVMCSample')])
nupdate_path = int(params[param_list.index('NExUpdatePath')])
seed = int(params[param_list.index('RndSeed')])
nsplit = int(params[param_list.index('NSplitSize')])
nstore = int(params[param_list.index('NStore')])
nsrcg = int(params[param_list.index('NSRCG')])

# write modpara.def
modpara.write('--------------------\n')
modpara.write('Model_Parameters  '+repr(0)+'\n')
modpara.write('--------------------\n')
modpara.write('VMC_Cal_Parameters\n')
modpara.write('--------------------\n')
modpara.write('CDataFileHead  '+datafile_head+'\n')
modpara.write('CParaFileHead  '+parafile_head+'\n')
modpara.write('--------------------\n')
modpara.write('NVMCCalMode    '+repr(calc_mode)+'\n')
modpara.write('NLanczosMode   '+repr(lanczos_mode)+'\n')
modpara.write('--------------------\n')
modpara.write('NDataIdxStart  '+repr(nidx_start)+'\n')
modpara.write('NDataQtySmp    '+repr(nsamp)+'\n')
modpara.write('--------------------\n')
modpara.write('Nsite          '+repr(Ns)+'\n')
modpara.write('Ncond          '+repr(ne)+'\n')
modpara.write('2Sz            '+repr(Sz)+'\n')
modpara.write('NSPGaussLeg    '+repr(ngauss)+'\n')
modpara.write('NSPStot        '+repr(nspstot)+'\n')
modpara.write('NMPTrans       '+repr(nmptrans)+'\n')
modpara.write('NSROptItrStep  '+repr(nSRsteps)+'\n')
modpara.write('NSROptItrSmp   '+repr(nSRsamp)+'\n')
modpara.write('DSROptRedCut   '+repr(SRredcut)+'\n')
modpara.write('DSROptStaDel   '+repr(SRstadel)+'\n')
modpara.write('DSROptStepDt   '+repr(SRdt)+'\n')
modpara.write('NVMCWarmUp     '+repr(nwarmup)+'\n')
modpara.write('NVMCInterval   '+repr(ninterval)+'\n')
modpara.write('NVMCSample     '+repr(nVMCsamp)+'\n')
modpara.write('NExUpdatePath  '+repr(nupdate_path)+'\n')
modpara.write('RndSeed        '+repr(seed)+'\n')
modpara.write('NSplitSize     '+repr(nsplit)+'\n')
modpara.write('NStore         '+repr(nstore)+'\n')
modpara.write('NSRCG          '+repr(nsrcg)+'\n')

modpara.close()

namelist = open('namelist.def', 'w')

namelist.write('\t'+'ModPara'+'\t'+'modpara.def'+'\n')
namelist.write('\t'+'LocSpin'+'\t'+'locspn.def'+'\n')
namelist.write('\t'+'Trans'+'\t'+'trans.def'+'\n')
namelist.write('\t'+'CoulombIntra'+'\t'+'coulombintra.def'+'\n')
namelist.write('\t'+'Orbital'+'\t'+'orbitalidx.def'+'\n')
namelist.write('\t'+'TransSym'+'\t'+'qptransidx.def'+'\n')
if use_jastrow:
    namelist.write('\t'+'Jastrow'+'\t'+'jastrowidx.def'+'\n')
namelist.write('\t'+'Gutzwiller'+'\t'+'gutzwilleridx.def'+'\n')
if input_fij:
    namelist.write('\t'+'InOrbital'+'\t'+'orbital.def'+'\n')
if read_fij:
    namelist.write('\t'+'InOrbital'+'\t'+'orbital.def'+'\n')

namelist.close()

# Write files for dvmc
run_dvmc = int(params[param_list.index('run_dvmc')])
if run_dvmc:

    nh = int(params[param_list.index('nh')])
    exc_per_site = int(params[param_list.index('exc_per_site')])

    # Write excitations first write spectrumpara.def
    dr1_x_min = 0
    dr1_x_max = Lx-1
    dr1_y_min = 0
    dr1_y_max = Ly-1
    dr2_x_min = 0
    dr2_x_max = Lx-1
    dr2_y_min = 0
    dr2_y_max = Ly-1

    eta = float(params[param_list.index('eta')])
    w_min = float(params[param_list.index('w_min')])
    w_max = float(params[param_list.index('w_max')])
    w_min_data = float(params[param_list.index('w_min_data')])
    w_max_data = float(params[param_list.index('w_max_data')])
    Nw = int(params[param_list.index('Nw')])
    exc_choice = str(params[param_list.index('exc_choice')])
    kpath_str = str(params[param_list.index('kPath')])

    spectrum = open('spectrumpara.def', 'w')
    spectrum.write("# input for postprocessing dvmc (dvmc.out with NVMCCalMode == 3)\n")
    spectrum.write("# '#' are treated as comment)\n")
    spectrum.write('\n')
    spectrum.write("# number of neighbors to include for makeExcitation_from_hopping.py\n")
    spectrum.write('nh           '+repr(nh)+'\n')
    spectrum.write('\n')
    spectrum.write("# number of excitations per site\n")
    spectrum.write('exc_per_site '+repr(exc_per_site)+'\n')
    spectrum.write('\n')
    spectrum.write("# parameters to generate excitation.def file:\n")
    spectrum.write("# these parameters are no longer used for makeExcitation_from_hopping.py\n")
    spectrum.write('dr1_x       '+repr(dr1_x_min)+':'+repr(dr1_x_max)+'\n')
    spectrum.write('dr1_y       '+repr(dr1_y_min)+':'+repr(dr1_y_max)+'\n')
    spectrum.write('dr2_x       '+repr(dr2_x_min)+':'+repr(dr2_x_max)+'\n')
    spectrum.write('dr2_y       '+repr(dr2_y_min)+':'+repr(dr2_y_max)+'\n')
    spectrum.write('\n')
    spectrum.write("# parameters to generate A(k,w) parameters\n")
    spectrum.write('eta          '+repr(eta)+'               # imaginary frequency\n')
    spectrum.write('w_min        '+repr(w_min)+'             # minimum frequency plotted\n')
    spectrum.write('w_max        '+repr(w_max)+'             # maximum frequency plotted\n')
    spectrum.write('w_min_data   '+repr(w_min_data)+'        # minimum frequency calculated\n')
    spectrum.write('w_max_data   '+repr(w_max_data)+'        # maximum frequency calculated\n')
    spectrum.write('Nw           '+repr(Nw)+'                # calculate '+repr(Nw)+' frequencies\n')
    spectrum.write('exc_choice   '+exc_choice+'        # take all excitations\n')
    spectrum.write('kPath        '+kpath_str+'               # k points chosen in the path to plot\n')

    spectrum.close()

    namelist_G = open('namelist_G.def', 'w')

    namelist_G.write('\t'+'ModPara'+'\t'+'modpara_G.def'+'\n')
    namelist_G.write('\t'+'LocSpin'+'\t'+'locspn.def'+'\n')
    namelist_G.write('\t'+'Trans'+'\t'+'trans.def'+'\n')
    namelist_G.write('\t'+'CoulombIntra'+'\t'+'coulombintra.def'+'\n')
    namelist_G.write('\t'+'Orbital'+'\t'+'orbitalidx.def'+'\n')
    namelist_G.write('\t'+'TransSym'+'\t'+'qptransidx.def'+'\n')
    namelist_G.write('\t'+'Jastrow'+'\t'+'jastrowidx.def'+'\n')
    namelist_G.write('\t'+'Gutzwiller'+'\t'+'gutzwilleridx.def'+'\n')
    namelist_G.write('\t'+'Excitation'+'\t'+'excitation.def'+'\n')

    namelist_G.close()

    modpara_G = open('modpara_G.def', 'w')

    # NVMCCalMode = 3 for dvmc
    calc_mode = 3

    # Write modpara_G.def
    modpara_G.write('--------------------\n')
    modpara_G.write('Model_Parameters  '+repr(0)+'\n')
    modpara_G.write('--------------------\n')
    modpara_G.write('VMC_Cal_Parameters\n')
    modpara_G.write('--------------------\n')
    modpara_G.write('CDataFileHead  '+'zvo'+'\n')
    modpara_G.write('CParaFileHead  '+parafile_head+'\n')
    modpara_G.write('--------------------\n')
    modpara_G.write('NVMCCalMode    '+repr(calc_mode)+'\n')
    modpara_G.write('--------------------\n')
    modpara_G.write('NDataIdxStart  '+repr(nidx_start)+'\n')
    modpara_G.write('NDataQtySmp    '+repr(nsamp)+'\n')
    modpara_G.write('--------------------\n')
    modpara_G.write('Nsite          '+repr(Ns)+'\n')
    modpara_G.write('Ncond          '+repr(ne)+'\n')
    modpara_G.write('2Sz            '+repr(Sz)+'\n')
    modpara_G.write('NSPGaussLeg    '+repr(ngauss)+'\n')
    modpara_G.write('NSPStot        '+repr(nspstot)+'\n')
    modpara_G.write('NMPTrans       '+repr(nmptrans)+'\n')
    modpara_G.write('NSROptItrStep  '+repr(nSRsteps)+'\n')
    modpara_G.write('NSROptItrSmp   '+repr(nSRsamp)+'\n')
    modpara_G.write('DSROptRedCut   '+repr(SRredcut)+'\n')
    modpara_G.write('DSROptStaDel   '+repr(SRstadel)+'\n')
    modpara_G.write('DSROptStepDt   '+repr(SRdt)+'\n')
    modpara_G.write('NVMCWarmUp     '+repr(nwarmup)+'\n')
    modpara_G.write('NVMCInterval   '+repr(ninterval)+'\n')
    modpara_G.write('NVMCSample     '+repr(nVMCsamp)+'\n')
    modpara_G.write('NExUpdatePath  '+repr(nupdate_path)+'\n')
    modpara_G.write('RndSeed        '+repr(seed)+'\n')
    modpara_G.write('NSplitSize     '+repr(nsplit)+'\n')
    modpara_G.write('NStore         '+repr(nstore)+'\n')
    modpara_G.write('NSRCG          '+repr(nsrcg)+'\n')

    modpara_G.close()
