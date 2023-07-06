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

R1x = [0, 1, 2, 3, 4, 5]
R1y = [0]
R2x = [-6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6]
R2y = [0]

L = W = 1

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
        if (param_name == 'Nc'):
            Nc = int(param_val)
else:
    print('no input file found!')
    sys.exit()

nSites = Nc + Nb


def find_neighbor_site(r, dx, dy):
    """Docs
    """
    r_x = (r+dx) % L
    r_y = (r/L+dy) % W
    r_out = (r_x + L*r_y) % nSites
    return r_out


def WriteExcitation():
    """Docs
    """
    # BEWARE: number of excitation per lattice site must be the
    # same for every sites.
    spectrumpara = open('spectrumpara.def').read()
    for line in spectrumpara.split('\n'):
        if len(line) > 0:
            if line[0] != '#':
                term = line.split()
                if term[0][:] == 'dr1_x':
                    dr1_x = ReadRange(term[1])
                if term[0][:] == 'dr1_y':
                    dr1_y = ReadRange(term[1])
                if term[0][:] == 'dr2_x':
                    dr2_x = ReadRange(term[1])
                if term[0][:] == 'dr2_y':
                    dr2_y = ReadRange(term[1])

    s = ''
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
        for dr1x in dr1_x:
            for dr1y in dr1_y:
                for dr2x in dr2_x:
                    for dr2y in dr2_y:
                        # condition to prevent unallowed excitation
                        if ((dr2x % L != 0) or (dr2y % W != 0)):
                            ra = find_neighbor_site(r1, dr1x, dr1y)
                            rb = find_neighbor_site(r1, dr2x, dr2y)
                            # condition to prevent redundant excitation
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


WriteExcitation()
