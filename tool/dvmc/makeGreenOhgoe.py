#!/usr/bin/env python3
import sys


def position(r):
    x = (r[0]+L[0]) % L[0]
    y = (r[1]+L[1]) % L[1]
    return (x, y)


def indexToPosition(i):
    # i = x+ y*Lx
    x = i % L[0]
    y = i/L[0]
    return (x, y)


def indexToWaveVector(k):
    # *****periodic boundary condition*****
    n = {}
    # change the index k to the corresponding grid
    n[0] = k % L[0]
    n[1] = k/L[0]
    # shift them to the first Brillouin zone (L/2:L/2]^d (Note: Lx=Ly is assumed)
    for d in range(Dim):
        if n[d] > L[d]/2:
            while n[d] > L[d]/2:
                n[d] -= L[d]
        else:
            pass
    kx = 2.0*pi*float(n[0])/L[0]
    ky = 2.0*pi*float(n[1])/L[1]

    return (kx, ky)


def indexToMinusWaveVector(k):
    # *****periodic boundary condition*****
    n = {}
    # change the index k to the corresponding grid
    n[0] = -(k % L[0])
    n[1] = -(k/L[0])
    # shift them to the first Brillouin zone (L/2:L/2]^d (Note: Lx=Ly is assumed)
    for d in range(Dim):
        if n[d] <= -L[d]/2:
            while n[d] <= -L[d]/2:
                n[d] += L[d]
        elif n[d] > L[d]/2:
            while n[d] > L[d]/2:
                n[d] -= L[d]
        else:
            pass

    kx = 2.0*pi*float(n[0])/L[0]
    ky = 2.0*pi*float(n[1])/L[1]

    return (kx, ky)


def positionToIndex(r):
    if r[0] < -L[0]:
        print("Error: r[0]<-Lx")
        sys.exit()
    if r[1] < -L[1]:
        print("Error: r[1]<-Ly")
        sys.exit()
    x = (r[0]+L[0]) % L[0]
    y = (r[1]+L[1]) % L[1]
    return x+L[0]*y


def positionToIndex2(r):
    # version 2 (return the sign for the antiperiodic condition)
    sgn = 1
    rx = r[0]
    ry = r[1]

    while rx < 0:
        rx += L[0]
        sgn *= -1
    while rx >= L[0]:
        rx -= L[0]
        sgn *= -1

    while ry < 0:
        ry += L[1]
    while ry >= L[1]:
        ry -= L[1]

    return (rx+L[0]*ry, sgn)


def neighborIndex(i, dr):
    r = indexToPosition(i)
    x = r[0]+dr[0]
    y = r[1]+dr[1]
    return positionToIndex([x, y])


def direction(i, j):
    ri = indexToPosition(i)
    rj = indexToPosition(j)
    dx = (rj[0]-ri[0]+L[0]) % L[0]
    dy = (rj[1]-ri[1]+L[1]) % L[1]
    return (dx, dy)


def locgrnIdx(i, j, s):

    return (Nsite*i+j)+s*Nsite*Nsite


def sgnAP(i, dr):
    r = indexToPosition(i)
    x = r[0]+dr[0]
    if (x >= L[0] or x <= -1):
        return -1
    else:
        return 1


def indexToDeltaVector(idx):
    # for d-wave SC correlation function
    # return a vector r which gives non-zero f_d(r).
    # f_d(r) is the form factor for the d-wave SC.

    if idx == 0:
        rx = 1
        ry = 0
    elif idx == 1:
        rx = -1
        ry = 0
    elif idx == 2:
        rx = 0
        ry = 1
    elif idx == 3:
        rx = 0
        ry = -1
    else:
        print("In deltaVector, i should be 0, 1, 2, or 3.")
        sys.exit()

    return (rx, ry)


def indexToFormFactor(idx):
    # return the value of the form factor f_d(r) (= 1 or -1)
    if idx == 0 or idx == 1:
        val = 1
    elif idx == 2 or idx == 3:
        val = -1
    else:
        print("In parity, idx should be 0, 1, 2, or 3.")
        sys.exit()

    return val


def vectorSum(r1, r2):
    rx = r1[0]+r2[0]
    ry = r1[1]+r2[1]

    return (rx, ry)


def vectorSum2(r1, r2, r3):
    rx = r1[0]+r2[0]+r3[0]
    ry = r1[1]+r2[1]+r3[1]

    return (rx, ry)


L = {}
L[0] = int(sys.argv[1])
L[1] = int(sys.argv[2])

Dim = 2
Nsite = L[0]*L[1]
eps = 1.0e-15
NSplitSize = 1

separator = '===============================\n'

# cisajs.def
f = open('greenone.def', 'w')
f.write(separator +
        "NCisAjs\t" + str(Nsite*Nsite*2) + "\n" +
        separator +
        "idx_i_j_s\n" +
        separator)
for i in range(Nsite):
    for j in range(Nsite):
        f.write("{0}\t{1}\t{2}\t{3}\n".format(i, 0, j, 0))
for i in range(Nsite):
    for j in range(Nsite):
        f.write("{0}\t{1}\t{2}\t{3}\n".format(i, 1, j, 1))
f.close()


# cisajscktalt.def
f = open('greentwo.def', 'w')
f.write(separator +
        "NCisAjsCktAlt\t" + str(Nsite*Nsite*(4+1+2+4*4)) + "\n" +
        # "NCisAjsCktAlt\t"+ str(0) +"\n"+
        separator +
        "===== i_s1_j_s2_k_t1_l_t2\n" +
        separator)
for i in range(Nsite):
    for j in range(Nsite):
        for s in range(2):
            for t in range(2):
                f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(
                    i, s, i, s, j, t, j, t))
for i in range(Nsite):
    for j in range(Nsite):
        f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(
            i, 0, j, 0, j, 1, i, 1))
        f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(
            i, 1, j, 1, j, 0, i, 0))
for i in range(Nsite):
    for j in range(Nsite):
        f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(
            i, 0, j, 0, i, 1, j, 1))

for i in range(Nsite):
    i1 = neighborIndex(i, [1, 0])
    i2 = neighborIndex(i, [-1, 0])
    i3 = neighborIndex(i, [0, 1])
    i4 = neighborIndex(i, [0, -1])
    for j in range(Nsite):
        j1 = neighborIndex(j, [1, 0])
        j2 = neighborIndex(j, [-1, 0])
        j3 = neighborIndex(j, [0, 1])
        j4 = neighborIndex(j, [0, -1])
        for s in range(1):
            for t in range(1):
                t = 1-t
                f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(
                    i, s, j, s, i1, t, j1, t))
                f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(
                    i, s, j, s, i1, t, j2, t))
                f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(
                    i, s, j, s, i1, t, j3, t))
                f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(
                    i, s, j, s, i1, t, j4, t))
                f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(
                    i, s, j, s, i2, t, j1, t))
                f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(
                    i, s, j, s, i2, t, j2, t))
                f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(
                    i, s, j, s, i2, t, j3, t))
                f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(
                    i, s, j, s, i2, t, j4, t))
                f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(
                    i, s, j, s, i3, t, j1, t))
                f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(
                    i, s, j, s, i3, t, j2, t))
                f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(
                    i, s, j, s, i3, t, j3, t))
                f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(
                    i, s, j, s, i3, t, j4, t))
                f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(
                    i, s, j, s, i4, t, j1, t))
                f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(
                    i, s, j, s, i4, t, j2, t))
                f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(
                    i, s, j, s, i4, t, j3, t))
                f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(
                    i, s, j, s, i4, t, j4, t))

f.close()
