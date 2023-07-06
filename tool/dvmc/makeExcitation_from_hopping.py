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

L, W = 1, 1

if os.path.isfile('StdFace.def'):
    print("StdFace.def found")
    StdFace = open('StdFace.def').read().replace(' ', '')
    if (StdFace.find('L=') >= 0):
        L = int(re.compile('L=([0-9]*)').findall(StdFace)[0])
    if (StdFace.find('W=') >= 0):
        W = int(re.compile('W=([0-9]*)').findall(StdFace)[0])
elif os.path.isfile('params'):
    print("using params file")
    paramfile = open('params', 'r')
    for line in paramfile:
        param_name = line.split()[0]
        param_val = line.split()[1]
        if (param_name == 'L'):
            L = int(param_val)
        if (param_name == 'W'):
            W = int(param_val)
        if (param_name == 'Nb'):
            Nb = int(param_val)
        if (param_name == 'U'):
            U = float(param_val)
else:
    print('no input file found!')
    sys.exit()

# load hopping matrix
if os.path.isfile('trans.def'):
    hop_mat = np.loadtxt('trans.def', skiprows=5)
else:
    print('no trans.def file found!')
    sys.exit()

# Cluster sites
nSites = L*W

# Bath sites
nSites_plus_Bath = int(hop_mat[:, 0].max()+1)

spectrumparaFileName = 'spectrumpara.def'
spectrumpara = open(spectrumparaFileName).read()

exc_per_site = 0
for line in spectrumpara.split('\n'):
    if len(line) > 0:
        if line[0] != '#':
            term = line.split()
            if term[0] == 'nh':
                nh = int(term[1])
            if term[0] == 'exc_per_site':
                exc_per_site = int(term[1])


def site_to_coords(site):
    x = int(site % W)
    y = int(site/W)
    return x, y


coords = np.zeros((L*W, 2))
for site in range(L*W):
    coords[site, :] = site_to_coords(site)


def find_neighbor_site_from_hopping(r, nh, hop_mat):
    # Array of neighbor sites for a given site
    neighbor_list = []

    # Array of hopping strengths to reach given neighbor site
    hop_list = np.zeros(1, dtype=np.float32)

    # Array of tuples of site pair and hoppings to reach site j from site i
    # i.e. [[site_i1, site_j1], [hop_11, hop_12, ...], [site_i2, site_j2, [hop_21, hop_22, ...]], ...]
    neighbor_plus_hops = []
    for nhop in range(nh):
        if (nhop == 0):
            for i in range(int(hop_mat.shape[0]/2)):
                coordsi = coords[int(hop_mat[i, 0]), :]
                coordsj = coords[int(hop_mat[i, 2]), :]
                dist = np.sqrt(np.sum((coordsi-coordsj)**2))
                if (hop_mat[i, 0] == r and dist < 2):
                    neighbor_list.append(int(hop_mat[i, 2]))
                    hop_list = [hop_mat[i, 4]]
                    new_elem = [[r, neighbor_list[-1]], hop_list]
                    neighbor_plus_hops.append(new_elem)
        else:
            old_neighbor_list = neighbor_list.copy()
            for rp in old_neighbor_list:
                rp_ind = old_neighbor_list.index(int(rp))
                for i in range(int(hop_mat.shape[0]/2)):

                    # Only include nearest, next-nearest neighbors
                    coordsi = coords[int(hop_mat[i, 0]), :]
                    coordsj = coords[int(hop_mat[i, 2]), :]
                    dist = np.sqrt(np.sum((coordsi-coordsj)**2))
                    if (hop_mat[i, 0] == rp and hop_mat[i, 2] != r and dist < 2 and hop_mat[i, 2] not in neighbor_list):
                        neighbor_list.append(int(hop_mat[i, 2]))
                        hop_list = neighbor_plus_hops[int(
                            rp_ind)][1]+[hop_mat[i, 4]]
                        new_elem = [[r, neighbor_list[-1]], hop_list]
                        neighbor_plus_hops.append(new_elem)

    return neighbor_plus_hops


def WriteExcitation():
    """Docs
    """
    # BEWARE: number of excitation per lattice site must
    # be the same for every sites.
    s = ''

    neighbor_list = []
    for r_ex in range(nSites):
        new_neighbors = find_neighbor_site_from_hopping(r_ex, nh, hop_mat)
        neighbor_list.append(new_neighbors)
        if (r_ex == 0):
            nn_max = len(new_neighbors)
        elif (len(new_neighbors) < nn_max):
            nn_max = len(new_neighbors)

    import copy
    # sort according to hoppings
    # for the moment this is just according to the number of hops
    # necessary to reach a given site, but should be changed in
    # the future
    for ns in range(len(neighbor_list)):
        weights = []
        for nhops in range(len(neighbor_list[ns])):
            w = 1.0/len(neighbor_list[ns][nhops][1])
            weights.append(w)
        w_inds = np.argsort(np.array(weights))[::-1]
        tmp = copy.deepcopy(neighbor_list[ns])
        for nhops in range(len(neighbor_list[ns])):
            neighbor_list[ns][nhops][:] = tmp[w_inds[nhops]]

    NExcitation = 0
    for r1 in range(nSites):

        s0 = ''
        s1 = ''
        s3 = ''
        s4 = ''
        # last two do not matter for type 0
        s0 += '%5d %4d %3d %3d \n' % (0, r1, 0, 0)
        # last one does not matter for type 1
        s1 += '%5d %4d %3d %3d \n' % (1, r1, r1, 0)
        NN = 2
        pairList_s3 = []
        pairList_s4 = []

        if (exc_per_site != 0):
            exc_ind = 0
            break_flag = False
            # while(exc_ind<exc_per_site or exc_ind>nn_max):
            for elem_a in neighbor_list[r1][:nn_max][:]:
                ra = elem_a[0][1]
                for elem_b in neighbor_list[r1][:nn_max][:]:
                    rb = elem_b[0][1]
                    # condition to prevent redundant excitation
                    if ((r1, ra, rb) not in pairList_s4):
                        NN += 1
                        exc_ind += 1
                        # everything matters for type 4
                        s4 += '%5d %4d %3d %3d \n' % (4, r1, ra, rb)
                        pairList_s4.append((r1, ra, rb))
                        if (exc_ind >= exc_per_site):
                            break_flag = True
                            break
                        # condition to prevent redundant excitation
                        if ((r1, rb, ra) not in pairList_s3):
                            NN += 1
                            # everything matters for type 4
                            s3 += '%5d %4d %3d %3d \n' % (3, r1, ra, rb)
                            pairList_s3.append((r1, ra, rb))
                            exc_ind += 1
                    if (exc_ind >= exc_per_site):
                        break_flag = True
                        break
                if (break_flag):
                    break
            s += s0 + s1 + s4 + s3
        else:
            # for ra in neighbor_list[r1][:nn_max]:
            for elem_a in neighbor_list[r1][:][:]:
                ra = elem_a[0][1]
                for elem_b in neighbor_list[r1][:][:]:
                    rb = elem_b[0][1]
                    # Condition to prevent redundant excitation
                    if ((r1, ra, rb) not in pairList_s4):
                        NN += 1
                        # everything matters for type 4
                        s4 += '%5d %4d %3d %3d \n' % (4, r1, ra, rb)
                        pairList_s4.append((r1, ra, rb))
                        # condition to prevent redundant excitation
                        if ((r1, rb, ra) not in pairList_s3):
                            NN += 1
                            # everything matters for type 4
                            s3 += '%5d %4d %3d %3d \n' % (3, r1, ra, rb)
                            pairList_s3.append((r1, ra, rb))
            s += s0 + s1 + s4 + s3

        # BEWARE: number of excitation per lattice site must be
        # the same for every sites.
        if r1 == 0:
            NExcitation = NN
        else:
            assert (NN == NExcitation)

    f = open('excitation.def', 'w')
    f.write(
        "===============================\n" +
        "NExcitation       "+str(NN)+"\n" +
        "L                 "+str(L)+"\n" +
        "W                 "+str(W)+"\n" +
        "----t---ri--ra--rb-------------\n"
    )

    f.write(s)
    f.close()


WriteExcitation()
