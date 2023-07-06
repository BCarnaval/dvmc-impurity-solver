#!/usr/bin/env python3
# zlib license:

# Copyright (c) 2019 Maxime Charlebois

# This software is provided 'as-is', without any express or implied
# warranty. In no event will the authors be held liable for any damages
# arising from the use of this software.

# Permission is granted to anyone to use this software for any purpose,
# including commercial applications, and to alter it and redistribute it
# freely, subject to the following restrictions:

# 1. The origin of this software must not be misrepresented; you must not
#    claim that you wrote the original software. If you use this software
#    in a product, an acknowledgment in the product documentation would be
#    appreciated but is not required.
# 2. Altered source versions must be plainly marked as such, and must not be
#    misrepresented as being the original software.
# 3. This notice may not be removed or altered from any source distribution.

import os
import re
import sys
import time
import numpy as np
from scipy.linalg import eigh
from numpy import linalg as la
import matplotlib.pyplot as plt

full_path = os.path.realpath(__file__)
pythonPathCode, file1 = os.path.split(full_path)

U, Nb = 0., 0
if os.path.isfile('StdFace.def'):
    StdFace = open('StdFace.def').read().replace(' ', '')
    if (StdFace.find('U=') >= 0):
        U = float(re.compile('U=([0-9.]*)').findall(StdFace)[0])
elif os.path.isfile('params'):
    paramfile = open('params', 'r')
    for line in paramfile:
        param_name = line.split()[0]
        param_val = line.split()[1]
        if (param_name == 'U'):
            U = float(param_val)
        if (param_name == 'Nb'):
            Nb = int(param_val)
        if (param_name == 'Nc'):
            Nc = int(param_val)
else:
    print("no input file found!")
    sys.exit()

read_qcm_params = False
if os.path.isfile('QCM_params.def'):
    read_qcm_params = True
    QCM_params = open('QCM_params.def')
    params = []
    for line in QCM_params:
        param_name = line.split()[0]
        param_val = line.split()[1]
        params.append([param_name, param_val])

trans_invariant = False
add_noise = False
qz_decomp = False

outputDir = 'output/'
spectrumparaFileName = 'spectrumpara.def'
verbose_read = 1

sum_rule_max_ok = 1.01
sum_rule_min_ok = 0.96

tol = 1e-10
k_tol = 3
pct_filter = 0.9

if (len(sys.argv) >= 2):
    spectrumparaFileName = sys.argv[1]
if (len(sys.argv) >= 3):
    outputDir = sys.argv[2]+'/'
if (len(sys.argv) >= 4):
    tol = float(sys.argv[3])
if (len(sys.argv) >= 5):
    addtl_filter = int(sys.argv[4])
if (len(sys.argv) >= 6):
    pct_filter = float(sys.argv[5])
if (len(sys.argv) >= 7):
    k_tol = float(sys.argv[6])


def dvmc_spectrum(verbose=1) -> None:
    """Docs
    """
    zqp_opt_dat = open(outputDir+'zqp_opt.dat').read()
    Omega = float((zqp_opt_dat.split())[0])
    if (read_qcm_params):
        params.append(['GS_energy', Omega])

    # Defaults
    w_min_data = -15.0
    w_max_data = 15.0
    eta = 0.2
    Nw = 2000
    exc_choice = [0, 1]
    spectrumpara = open(spectrumparaFileName).read()

    for line in spectrumpara.split('\n'):
        if len(line) > 0:
            if line[0] != '#':
                term = line.split()
                if term[0] == 'w_min_data':
                    w_min_data = float(term[1])
                if term[0] == 'w_max_data':
                    w_max_data = float(term[1])
                if term[0] == 'eta':
                    eta = float(term[1])
                if term[0] == 'Nw':
                    Nw = int(term[1])
                if (term[0][:] == 'exc_choice'):
                    if term[1][0:3] == 'all':
                        excitation_def = open('excitation.def').read()
                        line = excitation_def.split('\n')[1]
                        terms = line.split()
                        assert (terms[0] == 'NExcitation')
                        NExcitation = int(terms[1])
                        exc_choice = range(NExcitation)
                    elif term[1][0:6] == 'range(':
                        exc_choice = range(int(term[1][6:-1]))
                    else:
                        exc_choice = ReadRange(term[1])
                    if exc_choice[0] != 0:
                        print('error: first excitation must ALWAYS be 0.')
                        exit()

    if (verbose):
        print('excitations chosen:')
        print(exc_choice)

        print('\nOmega=', Omega)
        print('U=', U)

    len(exc_choice)
    Nc + Nb
    if (read_qcm_params):
        params.append(['Nc', Nc])
    dw = (w_max_data-w_min_data)/(Nw-1)
    w_ = np.array(range(Nw))*dw + w_min_data

    stot = time.time()
    S_CA = np.load(outputDir+'S_CA.npy')
    S_AC = np.load(outputDir+'S_AC.npy')
    H_CA = np.load(outputDir+'H_CA.npy')
    H_AC = np.load(outputDir+'H_AC.npy')

    spectrum_hole = np.zeros([Nc, Nc, Nw])
    spectrum_elec = np.zeros([Nc, Nc, Nw])

    # Symmetrize
    H_AC = 0.5*(H_AC+np.transpose(np.conjugate(H_AC)))
    S_AC = 0.5*(S_AC+np.transpose(np.conjugate(S_AC)))
    H_CA = 0.5*(H_CA+np.transpose(np.conjugate(H_CA)))
    S_CA = 0.5*(S_CA+np.transpose(np.conjugate(S_CA)))

    SV_AC, U_SVD_AC = eigh(S_AC)
    SV_CA, U_SVD_CA = eigh(S_CA)

    SV_AC_sorted = np.sort(SV_AC)[::-1]
    SV_CA_sorted = np.sort(SV_CA)[::-1]
    for ii in range(len(SV_AC)):
        print(ii, SV_AC_sorted[ii], SV_CA_sorted[ii])

    print('----------------------')

    # SVD to determine cutoff
    uu_ac, sing_v_ac, vh_ac = np.linalg.svd(S_AC, hermitian=True)
    uu_ca, sing_v_ca, vh_ca = np.linalg.svd(S_CA, hermitian=True)

    sing_v_ac_min = sing_v_ac[0]/(10**k_tol)
    reduced_inds = np.argwhere(sing_v_ac > sing_v_ac_min)
    n_keep = int(reduced_inds[-1])+1
    sing_v_ac_r = np.array(sing_v_ac[:n_keep])

    svals_ac_r = np.diag(sing_v_ac_r)
    uu_ac_r = uu_ac[:, :n_keep]
    vh_ac_r = vh_ac[:n_keep, :]

    S_AC = (uu_ac_r.dot(svals_ac_r)).dot(vh_ac_r)

    sing_v_ca_min = sing_v_ca[0]/(10**k_tol)
    reduced_inds = np.argwhere(sing_v_ca > sing_v_ca_min)
    n_keep = int(reduced_inds[-1])+1
    sing_v_ca_r = np.array(sing_v_ca[:n_keep])

    svals_ca_r = np.diag(sing_v_ca_r)
    uu_ca_r = uu_ca[:, :n_keep]
    vh_ca_r = vh_ca[:n_keep, :]

    S_CA = (uu_ca_r.dot(svals_ca_r)).dot(vh_ca_r)

    cond_num = np.log10(sing_v_ac_r[0]/sing_v_ac_r[-1])
    print("condition number S_AC = ", cond_num)

    cond_num = np.log10(sing_v_ca_r[0]/sing_v_ca_r[-1])
    print("condition number S_CA = ", cond_num)

    s = time.time()
    # Eigenvalue decomposition of overlap matrix, S
    SV_AC, U_SVD_AC = eigh(S_AC)
    SV_CA, U_SVD_CA = eigh(S_CA)
    e = time.time()
    print("Diagonalization 1 : ", e-s)

    SV_AC_sorted = np.sort(SV_AC)[::-1]
    SV_CA_sorted = np.sort(SV_CA)[::-1]

    np.abs(SV_AC.min())
    np.abs(SV_CA.min())

    nc_AC = len(SV_AC)
    nc_CA = len(SV_CA)

    break_AC = False
    break_CA = False
    SV_AC_tmp = SV_AC[::-1]
    SV_CA_tmp = SV_CA[::-1]

    for i in range(len(SV_AC)):
        if (SV_AC_tmp[i] < tol and not break_AC):
            nc_AC = i-1
            break_AC = True
        if (SV_CA_tmp[i] < tol and not break_CA):
            nc_CA = i-1
            break_CA = True
        if (break_AC and break_CA):
            break

    print("number of states kept AC: ", nc_AC)
    print("number of states kept CA: ", nc_CA)

    if (addtl_filter):
        nc_AC = int(pct_filter*nc_AC)
        nc_CA = int(pct_filter*nc_CA)

        print("number of states kept AC (w/addtl filtering): ", nc_AC)
        print("number of states kept CA (w/addtl filtering): ", nc_CA)

    s = time.time()
    D_AC_sqrt = np.sqrt(SV_AC[(S_AC.shape[0]-nc_AC):])
    b_AC = la.inv(U_SVD_AC)[(S_AC.shape[0]-nc_AC):, :]

    c = (D_AC_sqrt*b_AC.T).T
    S_AC_sqrt = U_SVD_AC[:, (S_AC.shape[0]-nc_AC):].dot(c)

    D_CA_sqrt = np.sqrt(SV_CA[(S_CA.shape[0]-nc_CA):])
    b_CA = la.inv(U_SVD_CA)[(S_CA.shape[0]-nc_CA):, :]
    c = (D_CA_sqrt*b_CA.T).T
    S_CA_sqrt = U_SVD_CA[:, (S_CA.shape[0]-nc_CA):].dot(c)

    D_AC_sqrt_inv = np.reciprocal(np.sqrt(SV_AC[(S_AC.shape[0]-nc_AC):]))
    c = (D_AC_sqrt_inv*b_AC.T).T
    Sinv_AC_sqrt = U_SVD_AC[:, (S_AC.shape[0]-nc_AC):].dot(c)

    D_CA_sqrt_inv = np.reciprocal(np.sqrt(SV_CA[(S_CA.shape[0]-nc_CA):]))
    c = (D_CA_sqrt_inv*b_CA.T).T
    Sinv_CA_sqrt = U_SVD_CA[:, (S_CA.shape[0]-nc_CA):].dot(c)

    e = time.time()
    print("Forming S, Sinv : ", e-s)

    s = time.time()
    Sinv_sqrt_H_Sinv_sqrt_AC = Sinv_AC_sqrt.dot(H_AC).dot(Sinv_AC_sqrt)
    Sinv_sqrt_H_Sinv_sqrt_CA = Sinv_CA_sqrt.dot(H_CA).dot(Sinv_CA_sqrt)
    e = time.time()
    print("Sinv*H*Sinv : ", e-s)

    s = time.time()
    e_ac, u_ac = eigh(Sinv_sqrt_H_Sinv_sqrt_AC)
    e_ca, u_ca = eigh(Sinv_sqrt_H_Sinv_sqrt_CA)
    e = time.time()
    print("Diagonalization 2 : ", e-s)
    e_ac - Omega
    -e_ca + Omega

    trim_inds = np.where(np.abs(e_ac) < 1e-10)
    e_ac_r = np.delete(e_ac, trim_inds)
    u_ac_r = np.delete(u_ac, trim_inds, axis=1)
    e_ac2_r = e_ac_r - Omega

    trim_inds = np.where(np.abs(e_ca) < 1e-10)
    e_ca_r = np.delete(e_ca, trim_inds)
    u_ca_r = np.delete(u_ca, trim_inds, axis=1)
    e_ca2_r = -e_ca_r + Omega

    s = time.time()

    # Compute Sbar^(1/2)*U_M
    u_ac = S_AC_sqrt.dot(u_ac_r)
    u_ca = S_CA_sqrt.dot(u_ca_r)
    e = time.time()
    print("Sbar^(1/2)*U_M : ", e-s)

    u_ac2 = u_ac[:Nc, :]
    u_ca2 = u_ca[:Nc, :]

    # Compute G(w) = Sbar^(1/2) * U_M * ((w + i * eta +- Omega) + E)^-1 * U_M^-1 * Sbar^(1/2)
    print("computing G(w) for range of frequencies")
    s = time.time()
    for ii in range(len(w_)):
        z_ac = w_[ii] + 1.j*eta + Omega
        z_ca = w_[ii] + 1.j*eta - Omega

        a = np.reciprocal(z_ac-e_ac_r)
        c1 = (a*u_ac2).T
        a = np.reciprocal(z_ca+e_ca_r)
        c2 = (a*u_ca2).T

        G_ac = u_ac2.dot(c1)  # O(x3)
        G_ca = u_ca2.dot(c2)

        spectrum_elec[:, :, ii] = -G_ac[:, :].imag/(np.pi)
        spectrum_hole[:, :, ii] = -G_ca[:, :].imag/(np.pi)

    totalAij = spectrum_hole + spectrum_elec
    e = time.time()
    print("Calculating spectrum : ", e-s)

    sum = (w_[1]-w_[0])*np.sum(totalAij[0, 0, :])
    e = time.time()
    print(sum)
    print("Execution time : ", e-stot)

    # Get dos
    dos = np.zeros([3, Nw], dtype='float')
    for i in range(Nw):
        dos[0, i] = np.trace(totalAij[:, :, i])
        dos[1, i] = np.trace(spectrum_elec[:, :, i])
        dos[2, i] = np.trace(spectrum_hole[:, :, i])
    file_dos = open(outputDir+'dos.dat', 'w')
    for ii in range(Nw):
        file_dos.write('% 7.6f   ' % w_[ii])
        for kk in range(3):
            file_dos.write('% 7.6f ' % (dos[kk, ii]))  # /(total_sum)))
        file_dos.write('\n')

    # Print solution for interface with QCM
    if (read_qcm_params):
        # Print_solution_for_QCM(params,u_ac2,e_ac2,u_ca2,e_ca2)
        print_solution_for_QCM(params, u_ac2, e_ac2_r, u_ca2, e_ca2_r)

    fig, ax = plt.subplots()

    shift = 0.5
    for ii in range(2):
        ax.plot(w_, totalAij[ii, ii, :]+ii*shift)
        ax.plot([w_.min(), w_.max()], [ii*shift, ii*shift], c='black')

    ax.axvline(x=0)
    ax.set_xlim(-6, 6)

    plt.savefig('spectrum_rspace.pdf')
    sys.exit()

    return


###############################################################################
# Sub routines
###############################################################################

def print_solution_for_QCM(params, u_ac2, e_ac2, u_ca2, e_ca2) -> None:
    """Docs
    """
    f = open(outputDir+'qmatrix.def', 'w')
    f2 = open(outputDir+'qmatrix.dat', 'w')

    for param in params:
        if (param[0] != 'sector' and param[0] != 'GS_energy' and param[0] != 'Nc'):
            f.write('%4s \t %4s\n' % (param[0].ljust(5), param[1]))
        if (param[0] == 'sector'):
            sector = param[1]
        if (param[0] == 'GS_energy'):
            GS_energy = param[1]
        if (param[0] == 'Nc'):
            Nc = param[1]

    nac = u_ac2.T.shape[0]
    nca = u_ca2.T.shape[0]

    Np = nac+nca

    f.write("\n")
    f.write('%10s  %3.8f %10s %10s\n' %
            ('GS_energy:', GS_energy, 'GS_sector:', sector+':1'))
    f.write('%10s %5s\n' % ('GF_format: ', 'bl'))
    f.write('%6s %10s\n' % ('mixing', '0'))
    f.write('state\n')
    f.write('%8s \t %3.8f \t %10s\n' % (sector, GS_energy, '1'))
    f.write('%2s \t %d \t %d\n' % ('w'.ljust(5), int(Nc), int(Np)))

    sortinds = np.argsort(e_ac2)
    e_ac2 = e_ac2[sortinds]
    u_ac2 = u_ac2[:, sortinds]

    sortinds = np.argsort(e_ca2)
    e_ca2 = e_ca2[sortinds]
    u_ca2 = u_ca2[:, sortinds]

    qmat = np.zeros((nac+nca, u_ac2.T.shape[1]+1))
    qmat[:nac, 0] = e_ac2
    qmat[:nac, 1:] = u_ac2.T
    qmat[nac:, 0] = e_ca2
    qmat[nac:, 1:] = u_ca2.T

    rows, cols = qmat.shape
    for i in range(rows):
        for j in range(cols):
            f.write("% 20.12e " % qmat[i, j])
            f2.write("% 20.12e " % qmat[i, j])
        f.write("\n")
        f2.write("\n")

    f.close()
    f2.close()

    return


def dotdot(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> np.ndarray:
    """Docs
    """
    return np.dot(np.dot(a, b), c)


def ReadRange(inputStr: str) -> list:
    rangeOut = []
    range1 = inputStr.split(',')
    for info in range1:
        element = info.split(':')
        if len(element) == 1:
            rangeOut += [int(element[0])]
        elif len(element) == 2:
            rangeOut += range(int(element[0]), int(element[1])+1)
        else:
            print('Your definition cannot be interpreted properly.')
            print('Please input numbers or pairs of numbers')
            print('sepated by comma (no spaces).')
            print('Pairs must contains only one ":"')
            print('example:')
            print('0:3,5,6,8:11')
            print("")
            print('will be interpreted as the list')
            print('[0,1,2,3,5,6,8,9,10,11]')
            exit()

    return rangeOut


if __name__ == "__main__":
    dvmc_spectrum(verbose_read)
