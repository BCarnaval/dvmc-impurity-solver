#!/usr/bin/env python3
import numpy as np


def sym():
    """Docs
    """
    flag_symmetrize_index = False  # f_ij  !=  f_ji
    symmetries = np.load('syms.npy')
    symmetries -= 1

    verbose = 1
    np.set_printoptions(precision=5)
    new_fij = find_f_ij(symmetries, flag_symmetrize_index, verbose)

    return new_fij

    exit()


###############################################################################
# End of the program
###############################################################################

def charHexa(digit: int) -> str:
    # Any cluster should not have above 64 sites, even in 2050 xD
    characters = '0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ!?'
    if digit < 64 and digit >= 0:
        char = characters[digit]
    else:
        char = '.'
    return char


def find_f_ij(symmetries, flag_symmetrize_index, verbose=2):
    gL = len(symmetries[0])
    print(gL)

    symmetryRuleList = symmetries

    def swapCouple(couple):
        return (couple[1], couple[0])

    f_ij_indexBasic = {}
    f_ij_index = {}
    for ii in range(0, gL):
        for jj in range(0, gL):
            f_ij_indexBasic[ii, jj] = ij(jj, ii)

            i = ii
            j = jj
            if flag_symmetrize_index:
                if ii > jj:
                    i = jj
                    j = ii

            f_ij_index[ii, jj] = ij(i, j)

    if verbose > 0:
        print('\ngeneral \nf_ij =')
        print(print_f_ij_indexToString(f_ij_indexBasic, gL))

        print('\ni->j, j->i symmetry: \nf_ij =')
        print(print_f_ij_indexToString(f_ij_index, gL))

    changed = True
    while changed:
        [changed, symmetryRuleList2] = symmetrizeElementOneByOne(
            symmetryRuleList, f_ij_index, gL, verbose)

    if verbose > 0:
        print('\nsymmetrized\nf_ij =')
        s, new_fij = print_f_ij_indexToString(
            f_ij_index, gL, return_fij_ind_array=True)
        print(s)

    fij_indep = np.zeros((gL, gL), dtype=int)
    listOfIndepGreenFound = {}
    N = 0
    for ii in range(0, gL):
        for jj in range(0, gL):
            if f_ij_index[ii, jj].string() not in listOfIndepGreenFound:
                N += 1
                listOfIndepGreenFound[f_ij_index[ii, jj].string()] = (ii, jj)
                fij_indep[ii, jj] = N-1
            else:
                fij_indep[ii, jj] = fij_indep[listOfIndepGreenFound[f_ij_index[ii, jj].string()]]

    print(N)
    print(fij_indep)

    return fij_indep, N


def symmetrizeElementOneByOne(symmetryRuleList, f_ij_index, gL, verbose):
    """Docs
    """
    ij(0, 0)
    symmetryRuleList2 = []

    for sym in symmetryRuleList:
        symmetryRuleList2.append(sym)
        symmetryRuleList2

        for i in range(0, gL):
            for j in range(0, gL):
                permutation_i = i
                permutation_j = j
                nn = 0
                isFirst = True

                # While didn't loop on orbit
                while ((permutation_i != i) or (permutation_j != j)) or isFirst:
                    nn += 1
                    permutation_i = sym[permutation_i]
                    permutation_j = sym[permutation_j]

                    if not (f_ij_index[permutation_j, permutation_i].isEqual(f_ij_index[j, i])):
                        ij0 = f_ij_index[j, i].string()
                        ij1 = f_ij_index[permutation_j,
                                         permutation_i].string()

                        # This portion of code ensure that always the lowest
                        # name (in alphabetical order) is kept
                        if ij0 < ij1:
                            f_ij_index[permutation_j, permutation_i].equal(
                                f_ij_index[j, i])
                        else:
                            f_ij_index[j, i].equal(
                                f_ij_index[permutation_j, permutation_i])

                        return [True, symmetryRuleList2]

                    isFirst = False
                    if nn > 100:
                        print('error, too many iteration')
                        exit()

        return [False, symmetryRuleList2]


class ij:
    def __init__(self, i, j):
        self.i = i
        self.j = j

    def equal(self, another):
        self.i = another.i
        self.j = another.j

    def isEqual(self, another):
        return (self.i == another.i) and (self.j == another.j)

    def string(self):
        return charHexa(self.i) + charHexa(self.j)


def print_f_ij_indexToString(f_ij_index, gL, return_fij_ind_array=False):
    """Docs
    """
    fij_new_ind = np.zeros((gL, gL, 2))
    maxLen = 1
    s1 = ''
    for jj in range(0, gL):
        for ii in range(0, gL):
            maxLen = max(maxLen, len(f_ij_index[ii, jj].string()))

    for jj in range(0, gL):
        for ii in range(0, gL):
            var = 1
            if (f_ij_index[ii, jj].string())[0] != '-':
                s1 += ' '
                var = 0
            s1 += f_ij_index[ii, jj].string()
            for kk in range(0, maxLen-len(f_ij_index[ii, jj].string())+var):
                s1 += ' '

        if jj is not gL-1:
            s1 += '\n'

    if (return_fij_ind_array):
        return s1, fij_new_ind
    else:
        return s1


if __name__ == "__main__":
    sym()
