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
import numpy as np
from ctypes import cdll, c_int, c_float, c_char_p

verbose = 0
full_path = os.path.realpath(__file__)
pythonPathCode, file1 = os.path.split(full_path)

if verbose:
    print(pythonPathCode)

Excitation = open('excitation.def').read().replace(' ', '')
n_exc = int(re.compile('NExcitation([0-9]*)').findall(Excitation)[0])

if verbose:
    print(n_exc)

Nc, Nb = 0, 0
if os.path.isfile('params'):
    paramfile = open('params', 'r')
    for line in paramfile:
        param_name = line.split()[0]
        param_val = line.split()[1]
        if (param_name == 'Nb'):
            Nb = int(param_val)
        if (param_name == 'Nc'):
            Nc = int(param_val)
else:
    print("no input file found!")
    sys.exit()


def main():
    dirOutput = './output/'

    fileIn = []
    n_file = len(sys.argv[:])-1
    n_file2 = 0

    if verbose:
        print(sys.argv[:])
        print(n_file)

    for nn in range(0, n_file):
        n_file2 += 1
        fileIn.append(sys.argv[1+nn])

    if verbose:
        print(fileIn)
    if (n_file < 1):
        print('error: no input files.\nexample\n$ mergeOutputBin.py output/zvo_nCHAm_nAHCm_00*\n')
        exit()

    Nsite = Nc + Nb
    print("Nsite = ", Nsite)
    is_so = True
    # use compiled version if present:
    if (not os.path.isfile(pythonPathCode+'/libdvmc_speedup.so')):
        if (not os.path.isfile(pythonPathCode+'/libdvmc_speedup.dylib')):
            print('Error: the shared object library "libdvmc_speedup.so" or "libdvmc_speedup.dylib" were not found.\nTo obtain it, go in the original directory of this python code (probably: ' + pythonPathCode+'/) and do make.')
            exit(0)
        is_so = False
    c_ArrayFloatN = c_float * (n_exc*n_exc*Nsite*Nsite)

    phys_CA_averaged = c_ArrayFloatN()
    phys_AC_averaged = c_ArrayFloatN()
    phys_CHA_averaged = c_ArrayFloatN()
    phys_AHC_averaged = c_ArrayFloatN()
    c_NExcitation = c_int(n_exc)
    c_dummy = c_int(0)
    c_Nsite = c_int(Nsite)
    c_int(n_exc*n_exc*Nsite*4)

    argList = (c_char_p * len(fileIn))()
    for i in range(len(fileIn)):
        fileIn[i] = fileIn[i].encode('utf-8')
    argList[:] = fileIn

    if verbose:
        print('***************')
        print(phys_CA_averaged[:10])

    # Loading and using our home made module:
    if is_so:
        lib1 = cdll.LoadLibrary(pythonPathCode+'/libdvmc_speedup.so')
    else:
        lib1 = cdll.LoadLibrary(pythonPathCode+'/libdvmc_speedup.dylib')
    lib1.mergeOutputBin(len(fileIn), argList, c_NExcitation, c_Nsite, c_dummy,
                        phys_CA_averaged, phys_AC_averaged, phys_CHA_averaged, phys_AHC_averaged, 0)

    def convert_c2numpy(phys_averaged):
        numpy_averaged = np.zeros(n_exc*n_exc*Nsite*Nsite)
        numpy_averaged[:] = phys_averaged[:]
        numpy_reshaped = numpy_averaged.reshape((Nsite*n_exc, Nsite*n_exc))
        return numpy_reshaped

    nCAm_up = convert_c2numpy(phys_CA_averaged)
    nACm_up = convert_c2numpy(phys_AC_averaged)
    nCHAm_up = convert_c2numpy(phys_CHA_averaged)
    nAHCm_up = convert_c2numpy(phys_AHC_averaged)

    np.set_printoptions(precision=2)
    print('\nS_AC')
    print(nACm_up)
    print('\nS_CA')
    print(nCAm_up)
    print('\nH_AC')
    print(nAHCm_up)
    print('\nH_CA')
    print(nCHAm_up)

    np.save(dirOutput+'S_CA', nCAm_up)
    np.save(dirOutput+'S_AC', nACm_up)
    np.save(dirOutput+'H_CA', nCHAm_up)
    np.save(dirOutput+'H_AC', nAHCm_up)


if __name__ == "__main__":
    main()
