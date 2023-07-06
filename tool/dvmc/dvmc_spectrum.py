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

import re
import os
import sys
import numpy as np
from ctypes import cdll
from scipy.linalg import eigh
from numpy import linalg as la
import matplotlib.pyplot as plt

full_path = os.path.realpath(__file__)
pythonPathCode, file1 = os.path.split(full_path)


StdFace = open('StdFace.def').read().replace(' ', '')
U, L, W = 0., 1, 1
if (StdFace.find('L=') >= 0):
    L = int(re.compile('L=([0-9]*)').findall(StdFace)[0])
if (StdFace.find('W=') >= 0):
    W = int(re.compile('W=([0-9]*)').findall(StdFace)[0])
if (StdFace.find('U=') >= 0):
    U = float(re.compile('U=([0-9.]*)').findall(StdFace)[0])

trans_invariant = False


def Xi(ri: int) -> int:
    """Docs
    """
    return ri % W


def Yi(ri: int) -> int:
    """Docs
    """
    return (ri//W) % L


outputDir = 'output/'
spectrumparaFileName = 'spectrumpara.def'
verbose_read = 1

sum_rule_max_ok = 1.01
sum_rule_min_ok = 0.96

if (len(sys.argv) >= 2):
    spectrumparaFileName = sys.argv[1]
if (len(sys.argv) >= 3):
    outputDir = sys.argv[2]+'/'
if (len(sys.argv) >= 4):
    verbose_read = int(sys.argv[3])
if (len(sys.argv) >= 5):
    sum_rule_min_ok = float(sys.argv[4])
if (len(sys.argv) >= 6):
    sum_rule_max_ok = float(sys.argv[5])
if (len(sys.argv) >= 7):
    print('example:\n$ vmc_spectrum.py \nor:\n$ vmc_spectrum.py spectrumpara'
          '.def\nor:\n$ vmc_spectrum.py spectrumpara.def output/')
    sys.exit()


def dvmc_spectrum(verbose=1):
    """Docs
    """
    zqp_opt_dat = open(outputDir+'zqp_opt.dat').read()
    Omega = float((zqp_opt_dat.split())[0])

    if (os.path.isfile(pythonPathCode+'/libdvmc_speedup.dylib')):
        cdll.LoadLibrary(pythonPathCode+'/libdvmc_speedup.dylib')
    else:
        print('Error: the shared object library "libdvmc_speedup.so" was not found.')
        print('To obtain it, go in the original directory of this python code')
        print('and compile the code with "make" (choose to link with LAPACK or MKL).')
        exit(-1)

    w_min_data = -15.0
    w_max_data = 15.0
    eta = 0.2
    Nw = 2000
    exc_choice = [0, 1]
    range(W * L)
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
                if term[0][:] == 'exc_choice':
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
        print('W=', W)
        print('L=', L)

    len(exc_choice)
    Nsite = W * L

    print("Nw: ", Nw)
    print("eta: ", eta)
    dw = (w_max_data-w_min_data)/(Nw-1)
    w_ = np.array(range(Nw))*dw + w_min_data

    S_CA = np.load(outputDir+'S_CA.npy')
    S_AC = np.load(outputDir+'S_AC.npy')
    H_CA = np.load(outputDir+'S_CHA.npy')
    H_AC = np.load(outputDir+'S_AHC.npy')

    np.zeros((Nw), np.cdouble)
    np.zeros((Nw), np.cdouble)

    spectrum_hole = np.zeros([Nsite, Nsite, Nw])
    spectrum_elec = np.zeros([Nsite, Nsite, Nw])

    if (not trans_invariant):

        # One big diagonalization [HS^-1] = Nsite*Nexc x Nsite*Nexc
        print("diagonalizing")
        e_ac, u_ac = eigh(H_AC, S_AC)
        e_ca, u_ca = eigh(H_CA, S_CA)

        print("computing U^-1")
        u_ac_m1 = la.inv(u_ac)
        u_ca_m1 = la.inv(u_ca)

        print("computing U^-1*S")
        us_ac = np.dot(u_ac_m1, S_AC)
        us_ca = np.dot(u_ca_m1, S_CA)

        len(e_ac)

        print("computing for each frequencies")
        for ii in range(len(w_)):
            z_ac = w_[ii] + 1.j*eta + Omega + U/2.
            z_ca = w_[ii] + 1.j*eta - Omega + U/2.

            G_ac = u_ac.dot(np.diag(np.reciprocal(z_ac-e_ac))
                            ).dot(us_ac)  # O(x3)
            G_ca = u_ca.dot(np.diag(np.reciprocal(z_ca+e_ca))).dot(us_ca)

            gac = G_ac[0:Nsite, 0:Nsite]
            gca = G_ca[0:Nsite, 0:Nsite]

            spectrum_elec[:, :, ii] = -gac[:, :].imag/(np.pi)
            spectrum_hole[:, :, ii] = -gca[:, :].imag/(np.pi)

        totalAij = spectrum_hole + spectrum_elec
        print("done computing A_ij(w)")

    fig, ax = plt.subplots()

    ax.plot(w_, totalAij[0, 0, :])
    ax.plot(w_, totalAij[1, 1, :]+0.0)
    ax.plot(w_, totalAij[2, 2, :]+0.0)
    ax.plot(w_, totalAij[3, 3, :]+0.0)
    ax.plot([-15, 15], [0, 0])

    plt.show()
    sys.exit()


###############################################################################
# Sub routines
###############################################################################

def dotdot(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> np.ndarray:
    """Docs
    """
    return np.dot(np.dot(a, b), c)


def ReadRange(inputStr: str) -> list:
    """Docs
    """
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
