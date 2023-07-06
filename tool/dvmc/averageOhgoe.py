#!/usr/bin/env python3

import sys
import glob
from math import *
import numpy as np


T1 = {
    (1, 0): 1.0,
    (-1, 0): 1.0,
    (0, 1): -1.0,
    (0, -1): -1.0
}
Tx = {
    (1, 0): 1.0,
    (-1, 0): 1.0
}
Ty = {
    (0, 1): 1.0,
    (0, -1): 1.0
}


def position(r):
    """Docs
    """
    x = (r[0] + Lx) % Lx
    y = (r[1] + Ly) % Ly
    return (x, y)


def InvVec(r):
    """Docs
    """
    return (-r[0], -r[1])


def indexToPosition(i):
    """Docs
    """
    x, y = i % Lx, i / Lx
    return (x, y)


def positionToIndex(r):
    """Docs
    """
    x = (r[0] + Lx) % Lx
    y = (r[1] + Ly) % Ly
    return x+Lx*y


def neighborIndex(i, dr):
    """Docs
    """
    r = indexToPosition(i)
    x = r[0] + dr[0]
    y = r[1] + dr[1]
    return positionToIndex([x, y])


def direction(i, j):
    """Docs
    """
    ri = indexToPosition(i)
    rj = indexToPosition(j)
    dx = (rj[0] - ri[0] + Lx) % Lx
    dy = (rj[1] - ri[1] + Ly) % Ly
    return (dx, dy)


def locgrnIdx(i, j, s):
    """Docs
    """
    return (Nsite * i + j) + s * Nsite * Nsite


def VecSum(ri, rj):
    """Docs
    """
    return (ri[0]+rj[0], ri[1]+rj[1])


def VecInv(ri):
    """Docs
    """
    return (-ri[0], -ri[1])


def VecDiff(ri, rj):
    """Docs
    """
    return (ri[0] - rj[0], ri[1] - rj[1])


def subIndex(i):
    """Docs
    """
    r = indexToPosition(i)
    sx = r[0] % Sx
    sy = r[1] % Sy
    return sx + (Sx * sy)


def subIndexToIndex(s):
    """Docs
    """
    sx = s % Sx
    sy = s / Sx
    i = positionToIndex([sx, sy])
    return i


def distance(r):
    """Returns the square of the distance between site i and site j.
    """
    rx, ry = r[0], r[1]
    if rx > Lx / 2:
        while rx > Lx / 2:
            rx -= Lx
    if ry > Ly / 2:
        while ry > Ly / 2:
            ry -= Ly

    return sqrt(rx * rx + ry * ry)


def sgnAP(i, dr):
    """Docs
    """
    if APFlag == 0:
        return 1.0

    else:
        r = indexToPosition(i)
        x = r[0] + dr[0]
        if (x < Lx - 0.5 and x > -0.5):
            return 1.0
        else:
            return -1.0


def WaveNumberList(Lx, Ly, APFlag):
    """Docs
    """
    wave = []
    if (APFlag == 1):
        for mx in xrange(-Lx / 2, Lx / 2):
            qx = 2.0 * pi * float(0.5 + mx) / float(Lx)
        for my in xrange(-Ly / 2, (Ly / 2) + 1):
            qy = 2.0*pi*float(my)/float(Ly)
            wave.append((qx, qy))
    else:
        for mx in range(-Lx / 2, (Lx / 2) + 1):
            qx = 2.0 * pi * float(mx) / float(Lx)
        for my in xrange(-Ly / 2, (Ly / 2) + 1):
            qy = 2.0 * pi * float(my) / float(Ly)
            wave.append((qx, qy))

    return wave


def CosSinList(Lx, Ly, APFlag):
    """Docs
    """
    cosList = {}
    sinList = {}
    Nsite = Lx*Ly
    for q in WaveNumberList(Lx, Ly, APFlag):
        (qx, qy) = q
        for ix in range(-Lx, Lx + 1):
            for iy in range(-Ly, Ly + 1):
                ri = (ix, iy)
                cosList[(q, ri)] = cos(qx * float(ix) + qy * float(iy))
                sinList[(q, ri)] = sin(qx * float(ix) + qy * float(iy))

    return [cosList, sinList]


def ReadPara(para):
    """Docs
    """
    ifile = open(para, 'r')
    for itr in xrange(5):
        ifile.readline().split()

    [name, DataFileHead] = ifile.readline().split()
    for itr in xrange(4):
        ifile.readline().split()

    [name, NDataIdxStart] = ifile.readline().split()
    [name, NDataQtySmp] = ifile.readline().split()
    ifile.readline().split()
    [name, Nsite] = ifile.readline().split()
    [name, Ncond] = ifile.readline().split()

    for itr in xrange(4):
        ifile.readline().split()

    [name, SRItrStep] = ifile.readline().split()
    [name, SRItrSmp] = ifile.readline().split()
    for itr in xrange(2):
        ifile.readline().split()

    [name, SRStepDt] = ifile.readline().split()
    ifile.close()

    Nsite = int(Nsite)
    Lx = int(sqrt(Nsite))
    Ly = int(sqrt(Nsite))
    Ne = int(Ncond) / 2
    NDataIdxStart = int(NDataIdxStart)
    NDataQtySmp = int(NDataQtySmp)
    SRItrStep = int(SRItrStep)
    SRItrSmp = int(SRItrSmp)
    SRStepDt = float(SRStepDt)

    params = [
        DataFileHead, NDataIdxStart, NDataQtySmp, Nsite, Ne, Lx, Ly,
        SRItrStep, SRItrSmp, SRStepDt
    ]

    return params


def GetFileList(flist):
    """Docs
    """
    ifile = open(flist, 'r')
    lines = []
    for line in ifile:
        line = line.rstrip("\n")
        lines.append(line)

    ifile.close()
    return lines


def ReadCisAjs(flist):
    """Docs
    """
    ifile = open(flist, 'r')
    ifile.readline().split()
    [name, NCisAjs] = ifile.readline().split()
    NCisAjs = int(NCisAjs)
    ciajList = []
    ciajListAppend = ciajList.append
    for itr in xrange(3):
        ifile.readline().split()

    for itr in xrange(NCisAjs):
        tmp = ifile.readline().split()
        tmp = map(int, tmp)
        [idx, i, j, s] = tmp
        ciajListAppend((i, j, s))

    ifile.close()

    return [NCisAjs, list(ciajList)]


def ReadCisAjsCktAlt(flist):
    """Docs
    """
    ifile = open(flist, 'r')
    ifile.readline().split()
    [name, NCisAjsCktAlt] = ifile.readline().split()
    NCisAjsCktAlt = int(NCisAjsCktAlt)
    ciajckalList = []
    ciajckalListAppend = ciajckalList.append
    for itr in xrange(3):
        ifile.readline().split()

    for itr in xrange(NCisAjsCktAlt):
        tmp = ifile.readline().split()
        tmp = map(int, tmp)
        [idx, idx2, i, j, s, k, l, t] = tmp
        ciajckalListAppend((i, j, s, k, l, t))

    ifile.close()

    return [NCisAjsCktAlt, list(ciajckalList)]


def ReadCisAjsCktAltDC(flist):
    """Docs
    """
    ifile = open(flist, 'r')
    ifile.readline().split()
    [name, NCisAjsCktAlt] = ifile.readline().split()
    NCisAjsCktAlt = int(NCisAjsCktAlt)
    ciajckalList = []
    ciajckalListAppend = ciajckalList.append
    for itr in xrange(3):
        ifile.readline().split()

    for itr in xrange(NCisAjsCktAlt):
        tmp = ifile.readline().split()
        tmp = map(int, tmp)
        [i, j, s, k, l, t] = tmp
        ciajckalListAppend((i, j, s, k, l, t))

    ifile.close()

    return [NCisAjsCktAlt, list(ciajckalList)]


###############################################################################
# Main : Calculate Physical Quantities
###############################################################################

if len(sys.argv) < 7:
    print("./average.py APFlag NkFlag Lx Ly Sx Sy")
    print("usage: APFlag==0 -> Periodic-Periodic boundary condition")
    print("       APFlag==1 -> AntiPeriodic-Periodic boundary condition\n")
    print("       NkFlag==0 -> Not Calculate momentum distribution n(k)")
    print("       NkFlag==1 -> Calculate momentum distribution n(k)\n")
    sys.exit()

APFlag = int(sys.argv[1])
NkFlag = int(sys.argv[2])
Lx = int(sys.argv[3])
Ly = int(sys.argv[4])
Sx = int(sys.argv[5])
Sy = int(sys.argv[6])

Nsub = Sx*Sy
Ntime = 1

# Read "*.def" files
paraList = ReadPara("../modpara_CA.def")

f = open('../greenone.def', 'r')
data = f.readline()
data = f.readline().split()
NCisAjs = int(data[1])
print(NCisAjs)
f.close()

f = open('../greentwo.def', 'r')
data = f.readline()
data = f.readline().split()
NCisAjsCktAlt = int(data[1])
print(NCisAjsCktAlt)
f.close()

pre = paraList[0]
Nsite = paraList[3]
Ne = paraList[4]
U = 1.0
pre2 = 'phy'
pre3 = 'ave'

WaveVector = WaveNumberList(Lx, Ly, APFlag)
[CosList, SinList] = CosSinList(Lx, Ly, APFlag)

if len(sys.argv) > 10:
    numList = sys.argv[10:]
else:
    numList = [n[len(pre) + 5:len(pre) + 8] for n in glob.glob(
        pre + "_out_*.dat")]

numList.sort()
nNum = len(numList)

for num in numList:
    # Read Green Functions
    ca = open(pre + '_cisajs_' + num + '.dat', 'r')
    caca = open(pre + '_cisajscktalt_' + num + '.dat', 'r')

    CisAjs = {}
    CisAjsCktAlt = {}
    time = [[] for i in xrange(Ntime)]

    # Read greenFunc
    for idx in xrange(NCisAjs):
        data = ca.readline().split()
        [i, s, j, t, val, val2] = data
        CisAjs[(int(i), int(j), int(s))] = float(val)

    # Read greenFunc2
    for idx in xrange(NCisAjsCktAlt):
        data = caca.readline().split()
        [i, s1, j, s2, k, t1, l, t2, val, val2] = data
        CisAjsCktAlt[(int(i), int(j), int(s1), int(k),
                      int(l), int(t1))] = float(val)

    ca.close()
    caca.close()

    # Calc Physical Quantities
    # Local Density
    ofile = open(pre2 + '_local_' + num + '.dat', 'w')
    ofile.write("# i  j  S_iS_j\n")
    for rx in xrange(Lx):
        for ry in xrange(Ly):
            i = rx + ry*Lx
            ofile.write("{0}\t{1}\t{2: .18f}\t{3: .18f}\n".format(
                rx, ry, CisAjs[(i, i, 0)].real, CisAjs[(i, i, 1)].real))

    ofile.close()

    # Spin Correlation
    print(num, "ss")
    ssdc = [[[[] for j in range(Nsite)] for i in xrange(Nsite)]
            for k in xrange(Ntime)]
    ss_xy = [[[[] for j in range(Nsite)]
              for i in xrange(Nsite)] for k in xrange(Ntime)]
    ss_z = [[[[] for j in range(Nsite)] for i in xrange(Nsite)]
            for k in xrange(Ntime)]
    ofile = open(pre2 + '_ss_' + num + '.dat', 'w')
    ofile.write("# i  j  S_iS_j\n")
    for k in xrange(Ntime):
        for i in xrange(Nsite):
            for j in xrange(Nsite):
                sz = 0.0
                sxy = 0.0
                if (i == j):
                    sz += 0.5*(CisAjs[(i, i, 0)] + CisAjs[(i, i, 1)])

                sxy += (-0.5) * (CisAjsCktAlt[(i, j, 0, j, i, 1)
                                              ] + CisAjsCktAlt[
                                                  (j, i, 0, i, j, 1)])
                sz += 0.25 * \
                    (CisAjsCktAlt[(i, i, 0, j, j, 0)] +
                     (CisAjsCktAlt[(i, i, 1, j, j, 1)]))
                sxy -= 0.25 * \
                    (CisAjsCktAlt[(i, i, 0, j, j, 1)] +
                     (CisAjsCktAlt[(i, i, 1, j, j, 0)]))
                s = (sxy + sz)
                ofile.write("{0}\t{1}\t{2: .18f}\t{3: .18f}\t{4: .18f}\t{5: .18f}\n".format(
                    i, j, s.real, s.imag, sz.real, sxy.real))
                ssdc[k][i][j] = s
                ss_xy[k][i][j] = sxy
                ss_z[k][i][j] = sz

    ofile.close()

    # Spin Structure Function
    print(num, "sq")
    ofile = open(pre2 + '_sq_' + num + '.dat', 'w')
    ofile.write(
        "# qx/pi qy/pi S(q).real  S(q).imag   Sz(q)/3  [Sx(q)+Sy(q)]/3 \n")
    for k in xrange(Ntime):
        mat = [[0.0 for ry in xrange(Ly)] for rx in xrange(Lx)]
        mat_xy = [[0.0 for ry in xrange(Ly)] for rx in xrange(Lx)]
        mat_z = [[0.0 for ry in xrange(Ly)] for rx in xrange(Lx)]

        for rx in xrange(Lx):
            for ry in xrange(Ly):
                for j in xrange(Nsite):
                    jx = j % Lx
                    jy = j / Lx
                    ix = (jx + rx) % Lx
                    iy = (jy + ry) % Ly
                    i = ix + (iy * Lx)
                    mat[rx][ry] += ssdc[k][i][j]
                    mat_xy[rx][ry] += ss_xy[k][i][j]
                    mat_z[rx][ry] += ss_z[k][i][j]

        sq = np.fft.fft2(mat)
        sq_xy = np.fft.fft2(mat_xy)
        sq_z = np.fft.fft2(mat_z)
        for mx in xrange(Lx + 1):
            qx = 2.0 * pi * float(mx) / float(Lx)
            for my in xrange(Ly + 1):
                qy = 2.0 * pi * float(my) / float(Ly)

                sqReal = sq[mx % Lx][my % Ly].real / (3.0 * float(Nsite))
                sqImag = sq[mx % Lx][my % Ly].imag / (3.0 * float(Nsite))
                sq_zReal = sq_z[mx % Lx][my % Ly].real / (3.0 * float(Nsite))
                sq_xyReal = sq_xy[mx % Lx][my % Ly].real / (3.0 * float(Nsite))
                ofile.write("{0}\t{1}\t{2: .18e}\t{3: .18e}\t{4: .18e}\t{5: .18e}\n".format(
                    qx / pi, qy / pi, sqReal, sqImag, sq_zReal, sq_xyReal))

            ofile.write("\n")

    ofile.close()

    # Momentum Distribution
    print(num, "nk")
    if (NkFlag == 1):
        ofile = open(pre2 + '_nk_' + num + '.dat', 'w')
        ofile.write("# kx/pi  ky/pi  n(k).real  n(k).imag\n")
        for k in xrange(Ntime):
            for q in WaveVector:
                (qx, qy) = q
                nkReal = 0.0
                nkImag = 0.0
                for i in xrange(Nsite):
                    ri = indexToPosition(i)
                    for j in xrange(Nsite):
                        rj = indexToPosition(j)
                        rij = VecDiff(ri, rj)

                        nkReal += CosList[(q, rij)] * \
                            (CisAjs[(i, j, 0)] + CisAjs[(i, j, 1)])
                        nkImag += SinList[(q, rij)] * \
                            (CisAjs[(i, j, 0)] + CisAjs[(i, j, 1)])

            nkReal /= 2.0*float(Nsite)
            nkImag /= 2.0*float(Nsite)
            ofile.write("{0: .10f}\t{1: .10f}\t{2: .18e}\t{3: .18e}\n".
                        format(qx / pi, qy / pi, nkReal.real, nkImag.imag))
            ofile.write("\n")

    ofile.close()

    # Density Correlation
    print(num, "nn")
    nndc = [[[[] for j in range(Nsite)] for i in xrange(Nsite)]
            for k in xrange(Ntime)]
    ofile = open(pre2 + '_nn_' + num + '.dat', 'w')
    ofile.write("# i  j  n_in_j\n")
    for k in xrange(Ntime):
        for i in xrange(Nsite):
            for j in xrange(Nsite):
                s = 0.0
                s += (CisAjsCktAlt[(i, i, 0, j, j, 0)] +
                      (CisAjsCktAlt[(i, i, 1, j, j, 1)]))
                s += (CisAjsCktAlt[(i, i, 0, j, j, 1)] +
                      (CisAjsCktAlt[(i, i, 1, j, j, 0)]))
                ofile.write("{0}\t{1}\t{2: .18f}\t{3: .18f}\n".format(
                    i, j, s.real, s.imag))
                nndc[k][i][j] = s

    ofile.close()

    # Charge Structure Factor
    print(num, "nq")
    ofile = open(pre2 + '_nq_' + num + '.dat', 'w')
    ofile.write("# qx qy  N(q).real  N(q).imag\n")
    for k in range(Ntime):
        mat = [[0.0 for ry in xrange(Ly)] for rx in xrange(Lx)]
        for rx in xrange(Lx):
            for ry in xrange(Ly):
                for j in xrange(Nsite):
                    jx = j % Lx
                    jy = j / Lx
                    ix = (jx + rx) % Lx
                    iy = (jy + ry) % Ly
                    i = ix + (iy * Lx)
                    mat[rx][ry] += CisAjsCktAlt[(i, i, 0, j, j, 0)] + \
                        CisAjsCktAlt[(i, i, 1, j, j, 1)]
                    mat[rx][ry] += CisAjsCktAlt[(i, i, 0, j, j, 1)] + \
                        CisAjsCktAlt[(i, i, 1, j, j, 0)]
                    mat[rx][ry] -= (2.0 * Ne / Nsite) * (2.0 * Ne / Nsite)

        nq = np.fft.fft2(mat)

        for mx in xrange(Lx + 1):
            qx = 2.0 * pi * float(mx) / float(Lx)
            for my in xrange(Ly + 1):
                qy = 2.0 * pi * float(my) / float(Ly)

                nqReal = nq[mx % Lx][my % Ly].real / float(Nsite)
                nqImag = nq[mx % Lx][my % Ly].imag / float(Nsite)
                ofile.write("{0: .10f}\t{1: .10f}\t{2: .18e}\t{3: .18e}\n".
                            format(qx / pi, qy / pi, nqReal, nqImag))
            ofile.write("\n")

    ofile.close()

    # Double Occupancy
    print(num, "do")
    ofile = open(pre2 + '_do_' + num + '.dat', 'w')
    ofile.write("# nn\n")

    for k in xrange(Ntime):
        nn = 0.0
        for isub in xrange(Nsub):
            i = subIndexToIndex(isub)
            nn += CisAjsCktAlt[(i, i, 0, i, i, 1)]
        nn /= float(Nsub)
        ofile.write("{0: .18e}\t{1: .18e}\n".format(nn.real, nn.imag))

    ofile.close()

    # Super Conducting Correlation Pd(ro)
    print(num, "dsc")
    ofile = open(pre2 + '_dsc_' + num + '.dat', 'w')
    ofile.write("# qx/pi qy/pi dSC(q).real  dSC(q).imag\n")
    ofile2 = open(pre2 + '_dsc2_' + num + '.dat', 'w')
    ofile2.write("# isub  rx   ry   P_dsc  \n")
    # Summed in terms of isub and (Lx-rx,Ly-ry)
    ofile3 = open(pre2 + '_dsc3_' + num+'.dat', 'w')
    ofile3.write("# rx   ry   P_dsc  \n")
    for t in xrange(Ntime):
        mat = [[0.0 for ry in xrange(Ly)] for rx in xrange(Lx)]
        SCmat = {}
        for rx in xrange(Lx):
            for ry in xrange(Ly):
                ro = [rx, ry]
                o = positionToIndex(ro)
                scAll = 0.0
                for i in xrange(Nsite):  # ri summation
                    sc = 0.0
                    io = neighborIndex(i, ro)
                    sgn_io = sgnAP(i, ro)
                    for dr_i, sgn_form_i in T1.iteritems():
                        i_dr = neighborIndex(i, dr_i)
                        sgn_i_dr = sgnAP(i, dr_i)
                        for dr_io, sgn_form_io in T1.iteritems():
                            io_dr = neighborIndex(io, dr_io)
                            ro_dr = [rx+dr_io[0], ry+dr_io[1]]
                            sgn_io_dr = sgnAP(i, ro_dr)
                            sgn_dsc = sgn_form_i * sgn_form_io
                            sgn_ap = sgn_io * sgn_i_dr * sgn_io_dr
                            sgn = sgn_dsc * sgn_ap

                            sc += CisAjsCktAlt[(i,    io,
                                                0, i_dr, io_dr, 1)] * sgn
                            sc += CisAjsCktAlt[(i,    io_dr,
                                                0, i_dr, io,    1)] * sgn
                            sc += CisAjsCktAlt[(i_dr, io,
                                                0, i,    io_dr, 1)] * sgn
                            sc += CisAjsCktAlt[(i_dr, io_dr,
                                                0, i,    io,    1)] * sgn

                            sc += CisAjsCktAlt[(io,    i,
                                                0, io_dr, i_dr, 1)] * sgn
                            sc += CisAjsCktAlt[(io,    i_dr,
                                                0, io_dr, i,    1)] * sgn
                            sc += CisAjsCktAlt[(io_dr, i,
                                                0, io,    i_dr, 1)] * sgn
                            sc += CisAjsCktAlt[(io_dr, i_dr,
                                                0, io,    i,    1)] * sgn

                            if i == io:
                                sc -= (CisAjs[(io_dr, i_dr, 0)] +
                                       CisAjs[(io_dr, i_dr, 1)]) * sgn
                                if i_dr == io_dr:
                                    sc += 2.0 * sgn
                            if i_dr == io_dr:
                                sc -= (CisAjs[(io, i, 0)] +
                                       CisAjs[(io, i, 1)]) * sgn
                            if i == io_dr:
                                sc -= (CisAjs[(io, i_dr, 0)] +
                                       CisAjs[(io, i_dr, 1)]) * sgn
                                if io == i_dr:
                                    sc += 2.0 * sgn
                            if io == i_dr:
                                sc -= (CisAjs[(io_dr, i, 0)] +
                                       CisAjs[(io_dr, i, 1)]) * sgn

                    SCmat[(rx, ry, i)] = sc / 4.0
                    scAll = scAll + sc
                mat[rx][ry] += scAll / (4.0 * Nsite)

        dscq = np.fft.fft2(mat)
        for mx in xrange(Lx + 1):
            qx = 2.0 * pi * float(mx) / float(Lx)
            for my in xrange(Ly + 1):
                qy = 2.0 * pi * float(my) / float(Ly)

                dscqReal = dscq[mx % Lx][my % Ly].real / float(Nsite)
                dscqImag = dscq[mx % Lx][my % Ly].imag / float(Nsite)
                ofile.write("{0: .10f}\t{1: .10f}\t{2: .18e}\t{3: .18e}\n".
                            format(qx / pi, qy / pi, dscqReal, dscqImag))
            ofile.write("\n")

        SCmat2 = [[[0.0 for isub in xrange(Nsub)]
                   for ry in xrange(Ly)] for rx in xrange(Lx)]
        for ix in xrange(Lx):
            for iy in xrange(Ly):
                i = positionToIndex((ix, iy))
                isub = subIndex(i)
                for rx in xrange(Lx):
                    for ry in xrange(Ly):
                        SCmat2[rx][ry][isub] += SCmat[(rx, ry, i)]

        for isub in range(Nsub):
            for rx in xrange(Lx):
                for ry in xrange(Ly):
                    ofile2.write(" {0}\t{1}\t{2}\t{3: .10e}\n".format(
                        isub, rx, ry, SCmat2[rx][ry][isub] / float(Nsite / Nsub)))

        for rx in xrange(Lx / 2 + 1):
            for ry in xrange(Ly / 2 + 1):
                ave = 0.0
                for isub in range(Nsub):
                    ave += SCmat2[rx][ry][isub]
                    rminus = indexToPosition(
                        positionToIndex((Lx - rx, Ly - ry)))
                    ave += SCmat2[rminus[0]][rminus[1]][isub]
                ofile3.write(" {0}\t{1}\t{2: .10e}\n".format(
                    rx, ry, ave / float(Nsite * 2)))

    ofile.close()
    ofile2.close()
    ofile3.close()

    # Super Conducting Correlation Ps(ro)
    print(num, "ssc")
    ofile = open(pre2 + '_ssc_' + num + '.dat', 'w')
    ofile.write("# i  j  Delta+_i_Delta_j\n")
    ofile2 = open(pre2 + '_ssc2_' + num + '.dat', 'w')
    ofile2.write("# i  j  Delta+_i_Delta_j\n")
    # Summed in terms of isub and (Lx-rx,Ly-ry)
    ofile3 = open(pre2 + '_ssc3_' + num + '.dat', 'w')
    ofile3.write("# rx   ry   P_ssc  \n")
    for t in xrange(Ntime):
        mat = [[0.0 for ry in xrange(Ly)] for rx in xrange(Lx)]
        SCmat = {}
        for rx in xrange(Lx):
            for ry in xrange(Ly):
                ro = [rx, ry]
                o = positionToIndex(ro)
                scAll = 0.0
                for i in xrange(Nsite):  # ri summation
                    sc = 0.0
                    io = neighborIndex(i, ro)
                    sc += CisAjsCktAlt[(i,  io,  0, i,  io, 1)]

                    if i == io:
                        sc += 1
                        sc -= (CisAjs[(io, i, 0)] + CisAjs[(io, i, 1)])

                    SCmat[(rx, ry, i)] = sc
                    scAll = scAll + sc

                mat[rx][ry] += scAll/(Nsite)

            maxP = {}
            for rx in xrange(Lx):
                for ry in xrange(Ly):
                    r = distance((rx, ry))
                    if r not in maxP:
                        maxP[r] = 0.0
                    if fabs(mat[rx][ry]) > fabs(maxP[r]):
                        maxP[r] = fabs(mat[rx][ry])

            numLine = 0
            for rx in xrange(Lx):
                for ry in xrange(Ly):
                    r = distance((rx, ry))
                    val = maxP[r]
                    if r <= sqrt(2.0) * float(Lx) * 0.5:
                        ofile.write("{0: .10e}\t{1: .10e}\t{2: .10e}\n".format(
                            r, val.real, val.imag))
                        numLine = numLine + 1

            SCmat2 = [
                [[0.0 for isub in xrange(Nsub)] for ry in xrange(Ly)] for rx in xrange(Lx)]

            for ix in xrange(Lx):
                for iy in xrange(Ly):
                    i = positionToIndex((ix, iy))
                    isub = subIndex(i)
                    for rx in xrange(Lx):
                        for ry in xrange(Ly):
                            SCmat2[rx][ry][isub] += SCmat[(rx, ry, i)]

            for isub in range(Nsub):
                for rx in xrange(Lx):
                    for ry in xrange(Ly):
                        ofile2.write(" {0}\t{1}\t{2}\t{3: .10e}\n".format(
                            isub, rx, ry, SCmat2[rx][ry][isub] / float(Nsite / Nsub)))

            for rx in xrange(Lx / 2 + 1):
                for ry in xrange(Ly / 2 + 1):
                    ave = 0.0
                    for isub in range(Nsub):
                        ave += SCmat2[rx][ry][isub]
                        rminus = indexToPosition(
                            positionToIndex((Lx - rx, Ly - ry)))
                        ave += SCmat2[rminus[0]][rminus[1]][isub]

                    ofile3.write(" {0}\t{1}\t{2: .10e}\n".format(
                        rx, ry, ave / float(Nsite * 2)))

    ofile.close()
    ofile2.close()
    ofile3.close()

###############################################################################
# Make Average Files
###############################################################################
ifile = [[] for i in xrange(nNum)]
ifile2 = [[] for i in xrange(nNum)]
# energy
print("ave", "ene")
for i in range(nNum):
    ifile[i] = open(pre + "_out_" + numList[i] + ".dat", 'r')

ofile = open(pre3 + "_ene.dat", 'w')
ofile.write("# E  sigma\n")
for k in xrange(Ntime):
    sum = [0.0, 0.0]
    for i in xrange(nNum):
        data = ifile[i].readline().split()
        time[k] = float(data[0])
        dtmp = float(data[0]) / float(Nsite)
        sum[0] += dtmp
        sum[1] += dtmp * dtmp

    ave = sum[0] / float(nNum)
    var = sum[1] / float(nNum) - ave * ave
    sigma = sqrt(abs(var / float(nNum)))

    ofile.write("{0: .10f}\t{1: .10f}\t{2: .10f}\n".format(
        time[k], ave, sigma))

for i in xrange(nNum):
    ifile[i].close()
ofile.close()

# local charge (spin) density
print("ave", "local")
for i in xrange(nNum):
    ifile[i] = open(pre2 + "_local_" + numList[i] + ".dat", 'r')
    data = ifile[i].readline()

ofile = open(pre3 + "_local.dat", 'w')
ofile.write("# i j S_iS_j sigma\n")
for k in xrange(Ntime):
    val = [[[] for i in xrange(Ly + 1)] for j in xrange(4)]
    for j in xrange(Nsite):
        sum = [0.0, 0.0]
        sum2 = [0.0, 0.0]
        for i in xrange(nNum):
            data = ifile[i].readline().split()
            itmp = int(data[0])
            jtmp = int(data[1])
            dtmp = float(data[2])
            dtmp2 = float(data[3])
            sum[0] += dtmp
            sum[1] += dtmp * dtmp
            sum2[0] += dtmp2
            sum2[1] += dtmp2 * dtmp2

        ave = sum[0] / float(nNum)
        var = sum[1] / float(nNum) - ave * ave
        ave2 = sum2[0] / float(nNum)
        var2 = sum2[1] / float(nNum) - ave2 * ave2
        sigma = sqrt(abs(var / float(nNum)))
        sigma2 = sqrt(abs(var2 / float(nNum)))
        ofile.write("{0}\t{1}\t{2: .10e}\t{3: .10e}\t{4: .10e}\t{5: .10e}\n".format(
            itmp, jtmp, ave, sigma, ave2, sigma2))
        if jtmp == Ly - 1:
            ofile.write("{0}\t{1}\t{2: .10e}\t{3: .10e}\t{4: .10e}\t{5: .10e}\n".format(
                itmp, jtmp + 1, ave, sigma, ave2, sigma2))
            ofile.write("\n")
        if itmp == Lx - 1:
            val[0][jtmp] = ave
            val[1][jtmp] = sigma
            val[2][jtmp] = ave2
            val[3][jtmp] = sigma2
            if jtmp == Ly - 1:
                val[0][jtmp + 1] = ave
                val[1][jtmp + 1] = sigma
                val[2][jtmp + 1] = ave2
                val[3][jtmp + 1] = sigma2
    if itmp == Lx - 1:
        for j in xrange(Ly + 1):
        ofile.write("{0}\t{1}\t{2: .10e}\t{3: .10e}\t{4: .10e}\t{5: .10e}\n".format(
            itmp + 1, j, val[0][j], val[1][j], val[2][j], val[3][j]))

    ofile.write("\n")

for i in xrange(nNum):
    ifile[i].close()
ofile.close()

# spin correlation ver.dc
print("ave", "ss")
for i in xrange(nNum):
    ifile[i] = open(pre2 + "_ss_" + numList[i] + ".dat", 'r')
    data = ifile[i].readline()

ofile = open(pre3 + "_ss.dat", 'w')
ofile.write("# i j S_iS_j sigma\n")
for k in xrange(Ntime):
    for j in xrange(Nsite * Nsite):
        sum = [0.0, 0.0]
        for i in xrange(nNum):
            data = ifile[i].readline().split()
            itmp = int(data[0])
            jtmp = int(data[1])
            dtmp = float(data[2])
            sum[0] += dtmp
            sum[1] += dtmp * dtmp

        ave = sum[0] / float(nNum)
        var = sum[1] / float(nNum) - ave * ave
        sigma = sqrt(abs(var/float(nNum)))
        ofile.write("{0}\t{1}\t{2: .10e}\t{3: .10e}\n".format(
            itmp, jtmp, ave, sigma))
        ofile.write("\n\n")

for i in xrange(nNum):
    ifile[i].close()
ofile.close()

# spin structure factor ver. dc
print("ave", "sq")
for i in xrange(nNum):
    ifile[i] = open(pre2 + "_sq_" + numList[i] + ".dat", 'r')
    data = ifile[i].readline()
ofile = open(pre3 + "_sq.dat", 'w')
ofile.write(
    "# qx/pi qy/pi  S(q).real sigma  S(q).imag sigma  Sz(q)/3 sigma [Sx(q)+Sy(q)]/3 sigma\n")
ofile2 = open(pre3 + "_sq_t.dat", 'w')
ofile2.write(
    "# time  S(pi,0).real sigma  S(pi,0).imag sigma S(pi,pi).real sigma S(pi,pi).imag sigma\n")
for k in xrange(Ntime):
    for mx in xrange(Lx + 1):
        for my in xrange(Ly + 1):
            sum = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            for i in xrange(nNum):
                data = ifile[i].readline().split()
                dtmpx = float(data[0])
                dtmpy = float(data[1])
                dtmp1 = float(data[2])
                dtmp2 = float(data[3])
                dtmp3 = float(data[4])
                dtmp4 = float(data[5])
                sum[0] += dtmp1
                sum[1] += dtmp1 * dtmp1
                sum[2] += dtmp2
                sum[3] += dtmp2 * dtmp2
                sum[4] += dtmp3
                sum[5] += dtmp3 * dtmp3
                sum[6] += dtmp4
                sum[7] += dtmp4 * dtmp4

            ave1 = sum[0] / float(nNum)
            var1 = sum[1] / float(nNum) - ave1 * ave1
            sigma1 = sqrt(abs(var1 / float(nNum)))
            ave2 = sum[2] / float(nNum)
            var2 = sum[3] / float(nNum) - ave2 * ave2
            sigma2 = sqrt(abs(var2 / float(nNum)))
            ave3 = sum[4] / float(nNum)
            var3 = sum[5] / float(nNum) - ave3 * ave3
            sigma3 = sqrt(abs(var3 / float(nNum)))
            ave4 = sum[6] / float(nNum)
            var4 = sum[7] / float(nNum) - ave4 * ave4
            sigma4 = sqrt(abs(var4 / float(nNum)))

            ofile.write("{0:.10f}\t{1:.10f}\t{2: .10e}\t{3: .10e}\t{4: .10e}\t{5: .10e}\t{6: .10e}\t{7: .10e}\t{8: .10e}\t{9: .10e}\n"
                        .format(dtmpx, dtmpy, ave1, sigma1, ave2, sigma2, ave3, sigma3, ave4, sigma4))
            if (dtmpx == 1.0) and (dtmpy == 0.0):
                ofile2.write("{0:.10e}\t{1: .10e}\t{2: .10e}\t{3: .10e}\t{4: .10e}\t"
                             .format(time[k], ave1, sigma1, ave2, sigma2, ave3, sigma3, ave4, sigma4))
            if (dtmpx == 1.0) and (dtmpy == 1.0):
                ofile2.write("{0:.10e}\t{1: .10e}\t{2: .10e}\t{3: .10e}\n"
                             .format(ave1, sigma1, ave2, sigma2, ave3, sigma3, ave4, sigma4))

        for i in xrange(nNum):
            ifile[i].readline()
        ofile.write("\n")

    ofile.write("\n\n")

for i in xrange(nNum):
    ifile[i].close()

ofile.close()
ofile2.close()

# Momentum distribution
if NkFlag == 1:
    print("ave", "nk")
    for i in range(nNum):
        ifile[i] = open(pre2 + "_nk_" + numList[i] + ".dat", 'r')
        data = ifile[i].readline()
    ofile = open(pre3 + "_nk.dat", 'w')
    ofile.write("# kx/pi ky/pi  n(k).real sigma  n(k).imag sigma\n")
    ofile2 = open(pre3 + "_dnk.dat", 'w')
    ofile2.write("# time  dnk(pi,0).real sigma\n")

    for k in xrange(Ntime):
        if (APFlag == 0):
            xList = xrange(-Lx / 2, Lx / 2 + 1)
        else:
            xList = xrange(-Lx / 2, Lx / 2)

        Ktotx = 0.0
        Ktoty = 0.0
        sigmax = 0.0
        sigmay = 0.0
        for mx in xList:
            for my in xrange(-Ly / 2, Ly / 2 + 1):
                sum = [0.0, 0.0, 0.0, 0.0]
                for i in xrange(nNum):
                    data = ifile[i].readline().split()
                    dtmpx = float(data[0])
                    dtmpy = float(data[1])
                    dtmp1 = float(data[2])
                    dtmp2 = float(data[3])
                    sum[0] += dtmp1
                    sum[1] += dtmp1 * dtmp1
                    sum[2] += dtmp2
                    sum[3] += dtmp2 * dtmp2

                ave1 = sum[0] / float(nNum)
                var1 = sum[1] / float(nNum) - ave1 * ave1
                sigma1 = sqrt(abs(var1 / float(nNum)))
                ave2 = sum[2] / float(nNum)
                var2 = sum[3] / float(nNum) - ave2 * ave2
                sigma2 = sqrt(abs(var2 / float(nNum)))
                Ktotx += dtmpx * ave1
                sigmax += dtmpx * dtmpx * sigma1 * sigma1
                Ktoty += dtmpy * ave1
                sigmay += dtmpy * dtmpy * sigma1 * sigma1
                ofile.write("{0:.10f}\t{1:.10f}\t{2: .10e}\t{3: .10e}\t{4: .10e}\t{5: .10e}\n"
                            .format(dtmpx, dtmpy, ave1, sigma1, ave2, sigma2))

                if mx == int(Lx / 2 - 1) and my == int(0):
                    ave3 = ave1
                    sigma3 = sigma1
                if mx == int(Lx / 2 - 1) and my == int(1):
                    ave4 = ave1
                    sigma4 = sigma1

            ofile.write("\n")
        ofile.write(" # Ktot_x  sigma  Ktot_y sigma \n")
        ofile.write(" \t{0: .10e}\t{1: .10e}\t{2: .10e}\t{3: .10e}\n"
                    .format(Ktotx, sqrt(sigmax), Ktoty, sqrt(sigmay)))
        for i in xrange(nNum):
            ifile[i].readline()

        ofile.write("\n\n")
        ofile2.write("{0:.10f}\t{1:.10f}\t{2: .10e}\n"
                     .format(time[k], ave3 - ave4, sqrt(sigma3 * sigma3 + sigma4 * sigma4)))

    for i in xrange(nNum):
        ifile[i].close()

    ofile.close()
    ofile2.close()


# Density correlation ver.dc
print("ave", "nn")
for i in xrange(nNum):
    ifile[i] = open(pre2 + "_nn_" + numList[i] + ".dat", 'r')
    data = ifile[i].readline()

ofile = open(pre3 + "_nn.dat", 'w')
ofile.write("# i j n_in_j sigma\n")
for k in xrange(Ntime):
    for i in xrange(Nsite):
        for j in xrange(Nsite):
            sum = [0.0, 0.0]
            for i in xrange(nNum):
                data = ifile[i].readline().split()
                itmp = int(data[0])
                jtmp = int(data[1])
                dtmp = float(data[2])
                sum[0] += dtmp
                sum[1] += dtmp * dtmp
                ave = sum[0] / float(nNum)
                var = sum[1] / float(nNum) - ave * ave
                sigma = sqrt(abs(var / float(nNum)))
                ofile.write("{0}\t{1}\t{2: .10e}\t{3: .10e}\n".format(
                    itmp, jtmp, ave, sigma))

        ofile.write("\n\n")

for i in xrange(nNum):
    ifile[i].close()
ofile.close()

# charge structure factor ver. dc
print("ave", "nq")
for i in xrange(nNum):
    ifile[i] = open(pre2 + "_nq_" + numList[i] + ".dat", 'r')
    data = ifile[i].readline()

ofile = open(pre3 + "_nq.dat", 'w')
ofile.write("# qx/pi qy/pi  N(q).real sigma  N(q).imag sigma\n")
ofile2 = open(pre3 + "_nq_t.dat", 'w')
ofile2.write(
    "# time  N(pi,0).real sigma  N(pi,0).imag sigma N(pi,pi).real sigma N(pi,pi).imag sigma\n")
for k in xrange(Ntime):
    for mx in xrange(Lx + 1):
        for my in xrange(Ly + 1):
            sum = [0.0, 0.0, 0.0, 0.0]
            for i in xrange(nNum):
                data = ifile[i].readline().split()
                dtmpx = float(data[0])
                dtmpy = float(data[1])
                dtmp1 = float(data[2])
                dtmp2 = float(data[3])
                sum[0] += dtmp1
                sum[1] += dtmp1 * dtmp1
                sum[2] += dtmp2
                sum[3] += dtmp2 * dtmp2

            ave1 = sum[0] / float(nNum)
            var1 = sum[1] / float(nNum) - ave1 * ave1
            sigma1 = sqrt(abs(var1 / float(nNum)))
            ave2 = sum[2] / float(nNum)
            var2 = sum[3] / float(nNum) - ave2 * ave2
            sigma2 = sqrt(abs(var2 / float(nNum)))
            ofile.write("{0:.10f}\t{1:.10f}\t{2: .10e}\t{3: .10e}\t{4: .10e}\t{5: .10e}\n"
                        .format(dtmpx, dtmpy, ave1, sigma1, ave2, sigma2))
            if (dtmpx == 1.0) and (dtmpy == 0.0):
                ofile2.write("{0:.10e}\t{1: .10e}\t{2: .10e}\t{3: .10e}\t{4: .10e}\t"
                             .format(time[k], ave1, sigma1, ave2, sigma2))

            if (dtmpx == 1.0) and (dtmpy == 1.0):
                ofile2.write("{0:.10e}\t{1: .10e}\t{2: .10e}\t{3: .10e}\n"
                             .format(ave1, sigma1, ave2, sigma2))

        for i in xrange(nNum):
            ifile[i].readline()

        ofile.write("\n")

    ofile.write("\n\n")

for i in xrange(nNum):
    ifile[i].close()
ofile.close()
ofile2.close()

# charge structure factor ver. dc
print("ave", "dscq")
for i in xrange(nNum):
    ifile[i] = open(pre2 + "_dsc_" + numList[i] + ".dat", 'r')
    data = ifile[i].readline()

ofile = open(pre3 + "_dscq.dat", 'w')
ofile.write("# qx/pi qy/pi dSC(q).real sigma  dSC(q).imag sigma\n")
for k in xrange(Ntime):
    for mx in xrange(Lx + 1):
        for my in xrange(Ly + 1):
            sum = [0.0, 0.0, 0.0, 0.0]
            for i in xrange(nNum):
                data = ifile[i].readline().split()
                dtmpx = float(data[0])
                dtmpy = float(data[1])
                dtmp1 = float(data[2])
                dtmp2 = float(data[3])
                sum[0] += dtmp1
                sum[1] += dtmp1 * dtmp1
                sum[2] += dtmp2
                sum[3] += dtmp2 * dtmp2

            ave1 = sum[0] / float(nNum)
            var1 = sum[1] / float(nNum) - ave1 * ave1
            sigma1 = sqrt(abs(var1 / float(nNum)))
            ave2 = sum[2] / float(nNum)
            var2 = sum[3] / float(nNum) - ave2 * ave2
            sigma2 = sqrt(abs(var2/float(nNum)))
            ofile.write("{0:.10f}\t{1:.10f}\t{2: .10e}\t{3: .10e}\t{4: .10e}\t{5: .10e}\n"
                        .format(dtmpx, dtmpy, ave1, sigma1, ave2, sigma2))

        for i in xrange(nNum):
            ifile[i].readline()
        ofile.write("\n")

for i in xrange(nNum):
    ifile[i].close()
ofile.close()

# double occupancy
print("ave", "do")
for i in xrange(nNum):
    ifile[i] = open(pre2 + "_do_" + numList[i] + ".dat", 'r')
    data = ifile[i].readline()
ofile = open(pre3 + "_do.dat", 'w')
ofile.write("# nn sigma\n")

for k in xrange(Ntime):
    sum = [0.0, 0.0]
    for i in xrange(nNum):
        data = ifile[i].readline().split()
        dtmp = float(data[0])
        sum[0] += dtmp
        sum[1] += dtmp * dtmp

    ave = sum[0] / float(nNum)
    var = sum[1] / float(nNum) - ave * ave
    sigma = sqrt(abs(var / float(nNum)))
    ofile.write("{0: .10f}\t{1: .10f}\t{2: .10f}\n".format(
        time[k], ave, sigma))

for i in xrange(nNum):
    ifile[i].close()
ofile.close()


# super conducting correlation ver.dc
print "ave", "dsc"
for i in xrange(nNum):
    ifile[i] = open(pre2 + "_dsc2_" + numList[i] + ".dat", 'r')
    data = ifile[i].readline()
for i in xrange(nNum):
    ifile2[i] = open(pre2 + "_dsc3_" + numList[i] + ".dat", 'r')
    data = ifile2[i].readline()

ofile = open(pre3 + "_dsc.dat", 'w')
ofile2 = open(pre3 + "_dsc2.dat", 'w')
# Correlation function along a particular direction
ofile3 = open(pre3 + "_dsc3.dat", 'w')
ofile.write("# r max|P_dsc(r)| sigma\n")
ofile2.write("# isub  rx  ry  P_dsc  sigma\n")
ofile3.write("# r   P(rx=r,0) sigma P(0,ry=r) sigma  r  P(rx,ry=rx) sigma\n")

for k in xrange(Ntime):
    phy = []
    print(numLine)

    # dsc2
    for j in xrange(Nsite * Nsub):
        sum = [0.0, 0.0]
        for n in xrange(nNum):
            data = ifile[n].readline().split()
            isub = int(data[0])
            rx = int(data[1])
            ry = int(data[2])
            tmp = float(data[3])
            sum[0] += tmp
            sum[1] += abs(tmp) * abs(tmp)

        ave = sum[0] / float(nNum)
        var = sum[1] / float(nNum) - abs(ave) * abs(ave)
        sigma = sqrt(abs(var / float(nNum)))
        list = [isub, rx, ry, ave, sigma]
        phy.append(list)

    for isub, rx, ry, ave, sigma in phy:
        ofile2.write(" {0}\t{1}\t{2}\t{3: .10e}\t{4: .10e}\n".format(
            isub, rx, ry, ave, sigma))

    # for dsc3
    ave_list = {}
    sigma_list = {}
    for j in xrange(((Lx / 2) + 1) * ((Ly / 2) + 1)):
        sum = [0.0, 0.0]
        for n in xrange(nNum):
            data = ifile2[n].readline().split()
            rx = int(data[0])
            ry = int(data[1])
            tmp = float(data[2])
            sum[0] += tmp
            sum[1] += abs(tmp) * abs(tmp)

        ave = sum[0] / float(nNum)
        var = sum[1] / float(nNum) - abs(ave) * abs(ave)
        sigma = sqrt(abs(var / float(nNum)))
        ave_list[(rx, ry)] = ave
        sigma_list[(rx, ry)] = sigma

    sc_ave = 0.0
    sc_sigma = 0.0
    icount = 0.0

    for rx in range(Lx/2+1):
        ofile3.write(" {0}\t{1: .10e}\t{2: .10e}\t{3: .10e}\t{4: .10e}\t{5: .10e}\t{6: .10e}\t{7: .10e}\n".format(rx, ave_list[(
            rx, 0)], sigma_list[(rx, 0)], ave_list[(0, rx)], sigma_list[(0, rx)], sqrt(2.0) * rx,  ave_list[(rx, rx)], sigma_list[(rx, rx)]))

    # for dsc
    maxP = {}
    err = {}
    for (rx, ry) in ave_list:
        r = distance((rx, ry))
        if r not in maxP:
            maxP[r] = 0.0
        if fabs(ave_list[(rx, ry)]) > fabs(maxP[r]):
            maxP[r] = fabs(ave_list[(rx, ry)])
            err[r] = sigma_list[(rx, ry)]

    sc_ave = 0.0
    sc_sigma = 0.0
    icount = 0.0
    rmax = (Lx / 2) * sqrt(2.0)

    for r in maxP:
        ave = maxP[r]
        sigma = err[r]
        if (r > 0) and (r < 10000):
            icount = icount + 1.0
            sc_ave = sc_ave + ave
            sc_sigma = sc_sigma + sigma * sigma
        ofile.write("{0: .10e}\t{1: .10e}\t{2: .10e}\n".format(r, ave, sigma))

    ofile.write("inf\t {0: .10e}\t{1: .10e}\n".format(
        sc_ave / icount, sqrt(sc_sigma) / icount))

for i in xrange(nNum):
    ifile[i].close()
for i in xrange(nNum):
    ifile2[i].close()

ofile.close()
ofile2.close()
ofile3.close()


# super conducting correlation ver.dc
print("ave", "ssc")
for i in xrange(nNum):
    ifile[i] = open(pre2 + "_ssc2_" + numList[i] + ".dat", 'r')
    data = ifile[i].readline()
for i in xrange(nNum):
    ifile2[i] = open(pre2 + "_ssc3_" + numList[i] + ".dat", 'r')
    data = ifile2[i].readline()

ofile = open(pre3+"_ssc.dat", 'w')
ofile2 = open(pre3+"_ssc2.dat", 'w')
# Correlation function along a particular direction
ofile3 = open(pre3+"_ssc3.dat", 'w')
ofile.write("# r max|P_ssc(r)| sigma\n")
ofile2.write("# isub  rx  ry  P_ssc  sigma\n")
ofile3.write("# r   P(rx=r,0) sigma P(0,ry=r) sigma  r  P(rx,ry=rx) sigma\n")
sum_tmp = 0.0
sigma2_tmp = 0.0

for k in xrange(Ntime):
    phy = []
    print(numLine)

    # ssc2
    for j in xrange(Nsite * Nsub):
        sum = [0.0, 0.0]
        for n in xrange(nNum):
            data = ifile[n].readline().split()
            isub = int(data[0])
            rx = int(data[1])
            ry = int(data[2])
            tmp = float(data[3])
            sum[0] += tmp
            sum[1] += abs(tmp) * abs(tmp)

        ave = sum[0] / float(nNum)
        var = sum[1] / float(nNum) - abs(ave) * abs(ave)
        sigma = sqrt(abs(var / float(nNum)))
        list = [isub, rx, ry, ave, sigma]
        phy.append(list)

    for isub, rx, ry, ave, sigma in phy:
        sum_tmp += ave
        sigma2_tmp += sigma * sigma
        ofile2.write(" {0}\t{1}\t{2}\t{3: .10e}\t{4: .10e}\n".format(
            isub, rx, ry, ave, sigma))

    ofile4 = open(pre3 + "_ps.dat", 'w')
    ofile4.write(
        "# P(0) sigma [ The definition of P(0) is Eq. (3.7) in Yokoyama, PTP 108, 59 (2002) ] \n")
    ofile4.write(" {0: .10e}\t{1: .10e}\n".format(
        sum_tmp/Nsub, sqrt(sigma2_tmp) / Nsub))
    ofile4.close()

    # ssc3
    ave_list = {}
    sigma_list = {}
    for j in xrange(((Lx / 2) + 1) * ((Ly / 2) + 1)):
        sum = [0.0, 0.0]
        for n in xrange(nNum):
            data = ifile2[n].readline().split()
            rx = int(data[0])
            ry = int(data[1])
            tmp = float(data[2])
            sum[0] += tmp
            sum[1] += abs(tmp) * abs(tmp)
        ave = sum[0] / float(nNum)
        var = sum[1] / float(nNum) - abs(ave) * abs(ave)
        sigma = sqrt(abs(var / float(nNum)))
        ave_list[(rx, ry)] = ave
        sigma_list[(rx, ry)] = sigma

    sc_ave = 0.0
    sc_sigma = 0.0
    icount = 0.0
    for rx in range((Lx / 2) + 1):
        ofile3.write(" {0}\t{1: .10e}\t{2: .10e}\t{3: .10e}\t{4: .10e}\t{5: .10e}\t{6: .10e}\t{7: .10e}\n".format(rx, ave_list[(
            rx, 0)], sigma_list[(rx, 0)], ave_list[(0, rx)], sigma_list[(0, rx)], sqrt(2.0) * rx,  ave_list[(rx, rx)], sigma_list[(rx, rx)]))

    # for ssc
    maxP = {}
    err = {}
    for (rx, ry) in ave_list:
        r = distance((rx, ry))
        if r not in maxP:
            maxP[r] = 0.0
        if fabs(ave_list[(rx, ry)]) > fabs(maxP[r]):
            maxP[r] = fabs(ave_list[(rx, ry)])
            err[r] = sigma_list[(rx, ry)]

    sc_ave = 0.0
    sc_sigma = 0.0
    icount = 0.0
    rmax = (Lx / 2) * sqrt(2.0)

    for r in maxP:
        ave = maxP[r]
        sigma = err[r]
        if (r > 0) and (r < 10000):
            icount = icount+1.0
            sc_ave = sc_ave + ave
            sc_sigma = sc_sigma + sigma*sigma
        ofile.write("{0: .10e}\t{1: .10e}\t{2: .10e}\n".format(r, ave, sigma))

    ofile.write("inf\t {0: .10e}\t{1: .10e}\n".format(
        sc_ave / icount, sqrt(sc_sigma) / icount))

for i in xrange(nNum):
    ifile[i].close()

for i in xrange(nNum):
    ifile2[i].close()

ofile.close()
ofile2.close()
ofile3.close()
