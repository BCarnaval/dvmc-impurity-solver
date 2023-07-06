#!/usr/bin/env python3
import numpy as np


def bound(i, i_max):
    """Docs
    """
    return i % i_max


def coor(L, dimen, index, n):
    """Docs
    """
    den = L
    coor = np.empty(dimen, dtype=int)
    for i in range(dimen-1, -1, -1):
        den /= n[i]
        coor[i] = index/den
        index %= den
    return coor


def index(coor, dimen, n):
    """Docs
    """
    lattice_index = 0
    den = 1
    for i in range(dimen):
        lattice_index += (coor[i]*den)
        den *= n[i]
    return lattice_index


def fold_vector(x, y, Wsub, Lsub):
    """Docs
    """
    fold_x = (x//Wsub)*Wsub
    fold_y = (y//Lsub)*Lsub

    return (fold_x, fold_y)


def fold(x, y, fold_x, fold_y, W, L):
    """Docs
    """
    new_x = (x - fold_x) % W
    new_y = (y - fold_y) % L
    new_index = new_x + new_y*W
    return new_index


# tb must be list with elements [i,j,t_ij]
# eb must be list with elements [i,i,eb_i]

def set_hopping(t, mu, tb, eb, Nc, Nb, Lx, Ly):
    """Docs
    """
    ndim = 2
    Ns = Nc + Nb
    L = [Lx, Ly]
    H0 = np.zeros((Ns, Ns), dtype=complex)

    # Set bath hopping
    if (Nb > 0):
        for i in range(len(tb[:, 0])):
            H0[int(tb[i, 0]), int(tb[i, 1])] = tb[i, 2]
            H0[int(tb[i, 1]), int(tb[i, 0])] = tb[i, 2]
        # Set bath on-site energy
        for i in range(len(eb[:, 0])):
            H0[int(eb[i, 0]), int(eb[i, 0])] = eb[i, 2]
    # Set cluster hopping and chemical potential
    for i in range(Nc):
        for j in range(Nc):
            # Add on-site chem potential or bath energy
            if (i == j):
                # Check if on cluster or bath
                H0[i, j] = mu
            else:
                coori = coor(Nc, ndim, i, L)
                coorj = coor(Nc, ndim, j, L)

                distx = bound(coorj[0]-coori[0], L[0])
                disty = bound(coorj[1]-coori[1], L[1])

                if (distx + disty == 1 or (distx == 0 and disty == Ly-1) or (disty == 0 and distx == Lx-1)):
                    H0[i, j] = t

    return H0

# Currently works only for the case of no spin-orbit coupling


def set_hopping_open(t, mu, tb, eb, Nc, Nb, Lx, Ly, tp=0.0, tpp=0.0):
    """Docs
    """
    ndim = 2
    Ns = Nc + Nb
    L = [Lx, Ly]
    H0 = np.zeros((Ns, Ns), dtype=complex)
    if (Nb > 0):
        for i in range(len(tb[:, 0])):
            H0[int(tb[i, 0]), int(tb[i, 1])] = tb[i, 2]
            H0[int(tb[i, 1]), int(tb[i, 0])] = tb[i, 2]
        # Set bath on-site energy
        for i in range(len(eb[:, 0])):
            H0[int(eb[i, 0]), int(eb[i, 0])] = eb[i, 2]
    # Set cluster hopping and chemical potential
    for i in range(Nc):
        for j in range(Nc):
            # Add on-site chem potential or bath energy
            if (i == j):
                # Check if on cluster or bath
                H0[i, j] = mu
            else:
                coori = coor(Nc, ndim, i, L)
                coorj = coor(Nc, ndim, j, L)

                dist = np.sum(np.abs(coori-coorj))
                distx = np.abs(coorj[0]-coori[0])
                disty = np.abs(coorj[1]-coori[1])

                if (np.abs(dist-1) < 1e-6):
                    H0[i, j] = t
                if (np.abs(distx-1) < 1e-6 and np.abs(disty-1) < 1e-6):
                    H0[i, j] = tp
                if ((np.abs(distx-2) < 1e-6 and np.abs(disty) < 1e-6) or (np.abs(disty-2) < 1e-6 and np.abs(distx) < 1e-6)):
                    H0[i, j] = tpp

    return H0


def fermi(energy, chem_poten, beta):
    """Docs
    """
    return np.exp(-beta*(energy-chem_poten))/(1+np.exp(-beta*(energy-chem_poten)))


def fij_from_H0(H0, ne, treat_open_shell, beta=30):
    evals, evecs = np.linalg.eigh(H0)

    for i in range(len(evals)):
        print(i+1, evals[i])

    degenerate = False
    # Check for gap b/w energies of ne and ne+1 levels
    # (i.e. whether or not they are degenerate)
    if (np.abs(evals[int(ne/2)-1]-evals[int(ne/2)]) < 1e-6):
        degenerate = True
        print("Open shell")

    ef = evals[int(ne/2)-1]
    print('Fermi Energy: ', ef)

    if (not treat_open_shell):
        degenerate = False

    fij = np.zeros((H0.shape[0]*H0.shape[1]), dtype=complex)
    ind = 0
    if (degenerate):
        for i in range(H0.shape[0]):
            for j in range(H0.shape[0]):
                for n in range(H0.shape[1]):
                    fij[ind] = fij[ind] + evecs[i, n] * \
                        evecs[j, n]*fermi(evals[n], ef, beta)
                ind += 1
    else:
        for i in range(H0.shape[0]):
            for j in range(H0.shape[0]):
                for n in range(int(ne/2)):
                    fij[ind] = fij[ind] + evecs[i, n]*evecs[j, n]
                ind += 1

    return fij


def fij_from_H_MF(H0, Lx, Ly, ne, U, avg_dens, delta_dens):
    """Docs
    """
    print("Getting fij from mean-field soln.")
    dimH = H0.shape[0]
    nup_MF = np.zeros(dimH)
    ndn_MF = np.zeros(dimH)

    for ix in range(Lx):
        for iy in range(Ly):
            nup_MF[ix+iy*Ly] = 0.5*U*(avg_dens-(-1)**(ix+iy)*delta_dens)
            ndn_MF[ix+iy*Ly] = 0.5*U*(avg_dens+(-1)**(ix+iy)*delta_dens)

    H0_MF = np.zeros((2*dimH, 2*dimH), dtype=complex)
    H0_MF[:dimH, :dimH] = H0 + np.diag(nup_MF)  # U*ndn_MF)
    H0_MF[dimH:, dimH:] = H0 + np.diag(ndn_MF)  # U*nup_MF)
    evals, evecs = np.linalg.eig(H0_MF)

    sort_inds = np.argsort(evals[:dimH])
    print(sort_inds)

    phi_up = evecs[:dimH, sort_inds]
    sort_inds = np.argsort(evals[dimH:])

    print(sort_inds+dimH)
    phi_dn = evecs[dimH:, sort_inds+dimH]

    fij = np.zeros((dimH**2), dtype=complex)
    ind = 0
    for i in range(dimH):
        for j in range(dimH):
            for n in range(int(ne/2)):
                fij[ind] = fij[ind] + phi_up[i, n]*phi_dn[j, n]
            ind += 1
    return fij


def write_hop_file(operator_list, model_file, Ns):
    """Docs
    """
    f = open(model_file, 'r')
    lines = f.readlines()

    op_vals = []
    for operator in operator_list:
        for index, line in enumerate(lines):
            ls = line.split()
            if (operator+' :' in line):
                op_val = float(line.split()[2])
                op_vals.append(op_val)

    hop_mat = []

    for operator in operator_list:
        op_ind = operator_list.index(operator)
        for index, line in enumerate(lines):
            ls = line.split()
            if (operator+' :' in line):
                op_val = float(line.split()[2])
            if (operator in ls and 'operator' in ls):
                lc = index+1
                while (True):
                    l = lines[lc]
                    lc += 1
                    l = l.split()
                    if (len(l)):
                        i = int(l[0])-1
                        j = int(l[1])-1
                        if (i < Ns and j < Ns):
                            if (len(hop_mat) == 0):
                                hop_mat = np.array(
                                    [i, 0, j, 0, op_vals[op_ind].real, op_vals[op_ind].imag])
                            else:
                                hop_mat = np.vstack(
                                    (hop_mat, np.array([i, 0, j, 0, op_vals[op_ind].real, op_vals[op_ind].imag])))
                    else:
                        break

    sort_inds = np.argsort(hop_mat[:, 0])
    hop_up = hop_mat[sort_inds]
    hop_mat[:, 1].fill(1)
    hop_mat[:, 3].fill(1)
    hop_dn = hop_mat[sort_inds]
    hop_mat = np.concatenate((hop_up, hop_dn))

    np.savetxt('hop.dat', hop_mat, fmt=[
               '%d', '%d', '%d', '%d', '%3.2f', '%3.2f'], delimiter='\t')

    return
