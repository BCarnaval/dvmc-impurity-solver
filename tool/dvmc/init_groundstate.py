#!/usr/bin/env python3
import sys
from math import (pi, cos, sin, sqrt)

###############################################################################
# Sub routines
###############################################################################


def position(r):
    x = (r[0]+L[0]) % L[0]
    y = (r[1]+L[1]) % L[1]
    return (x, y)


def indexToPosition(i):
    x = i % L[0]
    y = i/L[0]
    return (x, y)


def indexToWaveVectorBZ(k):
    n = {}
    # Change the index k to the corresponding grid
    n[0] = k % L[0]
    n[1] = k/L[0]

    # Move them to the first Brillouin zone (Note: Lx=Ly is assumed)
    for d in range(Dim):
        if n[d] <= -L[d]/2:
            while n[d] <= -L[d]/2:
                n[d] += L[d]
        elif n[d] > L[d]/2:
            while n[d] > L[d]/2:
                n[d] -= L[d]  # -= -Lx ## mistake found in 2014/9/4
        else:
            pass

    ky = 2.0*pi*float(n[1])/L[1]

    if APFlag == 1:
        # Antiperiodic boundary condition in the y-direction
        kx = 2.0*pi*float(n[0])/L[0] - (pi/L[0])
    else:
        kx = 2.0*pi*float(n[0])/L[0]

    return (kx, ky)


def waveVectorToIndexBZ(kwave):
    n = {}
    n[1] = round(kwave[1]*L[1]/(2.0*pi))
    if APFlag == 1:
        n[0] = round((kwave[0]+(pi/L[0]))*L[0]/(2.0*pi))
    else:
        n[0] = round(kwave[0]*L[0]/(2.0*pi))

    # Move them to the region [0:Lx-1][0:Ly-1]
    for d in range(Dim):
        if n[d] < 0:
            while n[d] < 0:
                n[d] += L[d]
        elif n[d] >= L[d]:
            while n[d] >= L[d]:
                n[d] -= L[d]
        else:
            pass

    return n[0]+n[1]*L[0]


def initKwaveBZ(kwaveBZ, kwaveIdxBZ):
    n = {}

    for k in range(Nsite):
        # Change the index k to the corresponding grid
        n[0] = k % L[0]
        n[1] = k/L[0]

        # Move them to the first Brillouin zone (Note: Lx=Ly is assumed)
        for d in range(Dim):
            if n[d] <= -L[d]/2:
                while n[d] <= -L[d]/2:
                    n[d] = n[d]+L[d]
            elif n[d] > L[d]/2:
                while n[d] > L[d]/2:
                    n[d] = n[d]-L[d]
                else:
                    pass

        if APFlag == 1:
            kwave = (2.0*pi*float(n[0])/L[0] -
                     (pi/L[0]), 2.0*pi*float(n[1])/L[1])
        else:
            kwave = (2.0*pi*float(n[0])/L[0], 2.0*pi*float(n[1])/L[1])
        kwaveBZ.append(kwave)
        kwaveIdxBZ[kwave] = k


def initKwaveAFBZ(kwaveAFBZ, kwaveIdxAFBZ):
    idx = 0
    n = {}

    for n[0] in range(L[0]):
        for n[1] in range(L[1]):

            # Move them to the first Brillouin zone
            for d in range(Dim):
                if n[d] <= -L[d]/2:
                    while n[d] <= -L[d]/2:
                        n[d] += L[d]
                elif n[d] > L[d]/2:
                    while n[d] > L[d]/2:
                        n[d] -= L[d]
                else:
                    pass

            # Check whether it belongs to the folded Brillouin zone
            kx_tmp = 2.0*pi*n[0]/float(L[0])
            ky_tmp = 2.0*pi*n[1]/float(L[1])

            if (kx_tmp+ky_tmp) <= (pi+eps) and (-kx_tmp+ky_tmp) <= (pi+eps) and (-kx_tmp+ky_tmp) > -pi+eps and (kx_tmp+ky_tmp) > -pi+eps:
                if APFlag == 1:
                    kwave = (kx_tmp - (pi/L[0]), ky_tmp)
                else:
                    kwave = (kx_tmp, ky_tmp)
                kwaveAFBZ.append(kwave)
                kwaveIdxAFBZ[kwave] = idx
                idx += 1
            else:
                pass

    if idx != (Nsite/2):
        print("Error in initKwaveAFBZ.\n")
        print("idx={0}\n".format(idx))
        sys.exit()


def positionToIndex(r):
    x = (r[0]+L[0]) % L[0]
    y = (r[1]+L[1]) % L[1]
    return x+L[0]*y


def direction(i, j):
    ri = indexToPosition(i)
    rj = indexToPosition(j)
    dx = (rj[0]-ri[0]+L[0]) % L[0]
    dy = (rj[1]-ri[1]+L[1]) % L[1]
    return (dx, dy)


def subIndex(i):
    r = indexToPosition(i)
    sx = r[0] % Sx
    sy = r[1] % Sy
    return sx+Sx*sy


def sgnAP(i, dr):
    r = indexToPosition(i)
    x = r[0]+dr[0]
    if (x >= L[0] or x <= -1):
        return -1
    else:
        return 1


###############################################################################
# Main calculations
###############################################################################
ifile = open(sys.argv[1], 'r')
ifile.readline().split()
ifile.readline().split()
ifile.readline().split()
ifile.readline().split()
ifile.readline().split()
[name, Lx] = ifile.readline().split()
[name, Ly] = ifile.readline().split()
[name, Sx] = ifile.readline().split()
[name, Sy] = ifile.readline().split()
[name, Ne] = ifile.readline().split()
[name, APFlag] = ifile.readline().split()
[name, CDWFlag] = ifile.readline().split()
[name, PhFlag] = ifile.readline().split()
[name, SymSc] = ifile.readline().split()
[name, GapSc] = ifile.readline().split()
[name, GapAF] = ifile.readline().split()
[name, T1] = ifile.readline().split()
[name, T2] = ifile.readline().split()
[name, T3] = ifile.readline().split()
[name, T4] = ifile.readline().split()
[name, T5] = ifile.readline().split()
ifile.close()

L = {}
L[0] = int(Lx)
L[1] = int(Ly)
Sx = int(Sx)
Sy = int(Sy)
Ne = int(Ne)
APFlag = int(APFlag)
CDWFlag = int(CDWFlag)
PhFlag = int(PhFlag)
SymSc = str(SymSc)
GapSc = float(GapSc)
GapAF = float(GapAF)
T1 = float(T1)
T2 = float(T2)
T3 = float(T3)
T4 = float(T4)
T5 = float(T5)

T2 = T2/T1
T3 = T3/T1
T4 = T4/T1
T5 = T5/T1
T1 = 1.0

Nsite = L[0]*L[1]
GapScA = {}
GapScB = {}
Dim = 2
Q = (pi, pi)
Disp = {}  # Energy dispersion
eps = 1.0e-13
Sgn = 1.0  # Change the sign if you want a CDW order
if CDWFlag == 1:
    Sgn = -1.0

kwaveBZ = []
kwaveIdxBZ = {}
kwaveAFBZ = []
kwaveIdxAFBZ = {}
initKwaveBZ(kwaveBZ, kwaveIdxBZ)
initKwaveAFBZ(kwaveAFBZ, kwaveIdxAFBZ)

# energy dispersion
for k in range(Nsite):  # BZ
    kwave = kwaveBZ[k]
    Disp[k] = -2.0 * T1 * (cos(kwave[0]) + cos(kwave[1])) \
        - 4.0 * T2 * cos(kwave[0])*cos(kwave[1]) \
        - 2.0 * T3 * (cos(2.0*kwave[0]) + cos(2.0*kwave[1])) \
        - 4.0 * T4 * (cos(2.0*kwave[0])*cos(kwave[1])
                      + cos(kwave[0])*cos(2.0*kwave[1])) \
        - 4.0 * T5 * cos(2.0*kwave[0])*cos(2.0*kwave[1])

if SymSc == "s":  # on-site s-wave
    for kAF in range(Nsite/2):  # AFBZ
        GapScA[kAF] = GapSc
        GapScB[kAF] = GapSc

elif SymSc == "2s":  # extended s-wave
    for kAF in range(Nsite/2):  # AFBZ
        kwave = kwaveAFBZ[kAF]
        GapScA[kAF] = GapSc * (cos(kwave[0]) + cos(kwave[1]))
        GapScB[kAF] = GapSc * (cos(kwave[0]+Q[0]) + cos(kwave[1]+Q[1]))

elif SymSc == "3s":  # extended s-wave
    for kAF in range(Nsite/2):  # AFBZ
        kwave = kwaveAFBZ[kAF]
        GapScA[kAF] = GapSc * cos(kwave[0]) * cos(kwave[1])
        GapScB[kAF] = GapSc * cos(kwave[0]+Q[0]) * cos(kwave[1]+Q[1])

elif SymSc == "2d":  # d_{x^2-y^2}-wave
    for kAF in range(Nsite/2):  # AFBZ
        kwave = kwaveAFBZ[kAF]

        # We should impose the relation GapScB(k)=GapScA(k+Q).
        # Otherwise, in the limit of GapAF->0, our wave function does not
        # reduce to the BCS wave function.
        GapScA[kAF] = GapSc * (cos(kwave[0]) - cos(kwave[1]))
        GapScB[kAF] = GapSc * (cos(kwave[0]+Q[0]) - cos(kwave[1]+Q[1]))
else:
    print("Temporarily, only s-wave or d-wave are available for the symmetry of superconductors.\n")
    sys.exit()

# energy dispersion of AF bands [Eq.(6)]
kwavePlusQ = {}
xi1 = {}
xi2 = {}
dispAF_a = {}
dispAF_b = {}
for kAF in range(Nsite/2):  # AFBZ
    kwave = kwaveAFBZ[kAF]
    kwavePlusQ[0] = kwave[0] + Q[0]
    kwavePlusQ[1] = kwave[1] + Q[1]
    kPlusQ = waveVectorToIndexBZ(kwavePlusQ)
    k = waveVectorToIndexBZ(kwave)

    xi1[kAF] = (Disp[k] - Disp[kPlusQ]) * 0.5
    xi2[kAF] = (Disp[k] + Disp[kPlusQ]) * 0.5
    dispAF_a[kAF] = xi2[kAF] - sqrt(xi1[kAF]**2 + GapAF**2)
    dispAF_b[kAF] = xi2[kAF] + sqrt(xi1[kAF]**2 + GapAF**2)

# determine the chemical potential
etmp = [[] for i in range(Nsite*Nsite)]
icount = 0
for kAF in range(Nsite/2):  # AFBZ
    etmp[icount] = sqrt(dispAF_a[kAF]**2+GapScA[kAF]**2)
    icount = icount+1
    etmp[icount] = -sqrt(dispAF_a[kAF]**2+GapScA[kAF]**2)
    icount = icount+1
    etmp[icount] = sqrt(dispAF_b[kAF]**2+GapScB[kAF]**2)
    icount = icount+1
    etmp[icount] = -sqrt(dispAF_b[kAF]**2+GapScB[kAF]**2)
    icount = icount+1
etmp.sort()


mu = (etmp[2*Ne-1]+etmp[2*Ne])*0.5
print("mu={0}".format(mu))

# BCS coefficients [Eq.(8)]
u = {}
v = {}
for kAF in range(Nsite/2):  # AFBZ
    denom = sqrt(xi1[kAF]**2+GapAF**2)
    if abs(denom) < eps:
        print("sqrt( xi1(k)^2+GapAF^2 ) is too small or zero.\n")
        sys.exit()
    tmp = xi1[kAF]/denom  # Caution: the denominator becomes zero if GapAF=0.0
    u[kAF] = sqrt(0.5 * (1-tmp))
    v[kAF] = sqrt(0.5 * (1+tmp))

# Coefficient in the k-representation of the wave function
# (in terms of bogoliubov particles) [Eq.(9)]
psiA = {}
psiB = {}
psi1 = {}
psi2 = {}
for kAF in range(Nsite/2):  # AFBZ
    tmp1 = dispAF_a[kAF] - mu
    tmp2 = dispAF_b[kAF] - mu
    denom1 = tmp1 + sqrt(tmp1**2 + GapScA[kAF]**2)
    denom2 = tmp2 + sqrt(tmp2**2 + GapScB[kAF]**2)

    psiA[kAF] = GapScA[kAF]/(denom1+1.0e-5)
    psiB[kAF] = GapScB[kAF]/(denom2+1.0e-5)

# Coefficient in the k-representation of the wave function [Eq.(16)]
for kAF in range(Nsite/2):  # AFBZ
    a = (u[kAF]**2)*psiA[kAF] - Sgn * (v[kAF]**2)*psiB[kAF]
    b = - Sgn * (v[kAF]**2)*psiA[kAF] + (u[kAF]**2)*psiB[kAF]
    c = (Sgn*psiA[kAF] + psiB[kAF]) * u[kAF] * v[kAF]

    kwave = kwaveAFBZ[kAF]
    kwavePlusQ[0] = kwave[0] + Q[0]
    kwavePlusQ[1] = kwave[1] + Q[1]
    k = waveVectorToIndexBZ(kwave)
    kPlusQ = int(waveVectorToIndexBZ(kwavePlusQ))
    psi1[k] = a
    psi1[kPlusQ] = b
    psi2[kAF] = c


# Coefficient in the r-representation of the wave function [Eq.(19)]
fRe = {}
fIm = {}
for i in range(Nsite):
    for j in range(Nsite):

        fRe[i+Nsite*j] = 0.0
        fIm[i+Nsite*j] = 0.0
        ri = indexToPosition(i)
        rj = indexToPosition(j)

        for k in range(Nsite):  # BZ
            kwave = kwaveBZ[k]
            phase = kwave[0]*(ri[0]-rj[0]) + kwave[1]*(ri[1]-rj[1])
            fRe[i+Nsite*j] += psi1[k]*cos(phase)
            fIm[i+Nsite*j] += psi1[k]*sin(phase)

        for kAF in range(Nsite/2):  # AFBZ
            kwave = kwaveAFBZ[kAF]
            phase1 = kwave[0]*(ri[0]-rj[0]) + kwave[1]*(ri[1]-rj[1])
            phase2 = phase1 + Q[0]*ri[0] + Q[1]*ri[1]
            phase3 = phase1 - (Q[0]*rj[0] + Q[1]*rj[1])
            fRe[i+Nsite*j] += psi2[kAF] * (cos(phase2) - Sgn * cos(phase3))
            fIm[i+Nsite*j] += psi2[kAF] * (sin(phase2) - Sgn * sin(phase3))

        fRe[i+Nsite*j] = fRe[i+Nsite*j]/Nsite
        fIm[i+Nsite*j] = fIm[i+Nsite*j]/Nsite

# orbitalIdx file
fname = 'orbitalidx.def'
Nidx = 0
Nsub = Sx*Sy
NOrbital = Nsub*Nsite
orbital = [[-1 for dj in range(Nsite)] for isub in range(Nsub)]
orbitalSiteI = {}
orbitalSiteJ = {}

f = open(fname, 'w')
separator = '---------------------------------------------\n'
f.write(separator +
        "NOrbitalIdx\t" + str(NOrbital) + "\n" +
        "ComplexType\t" + str(0) + "\n" +
        separator +
        separator)

for isub in range(Nsub):
    for dj in range(Nsite):
        orbital[isub][dj] = Nidx
        Nidx += 1

for i in range(Nsite):
    for j in range(Nsite):
        isub = subIndex(i)
        dj = positionToIndex(direction(i, j))
        idx = orbital[isub][dj]
        orbitalSiteI[idx] = i
        orbitalSiteJ[idx] = j
        f.write("{0}\t{1}\t{2}\n".format(i, j, orbital[isub][dj]))

for i in range(Nidx):
    f.write("{0}\t{1}\n".format(i, 1))
f.close

# Read files
gut = open('gutzwilleridx.def', 'r')
data = gut.readline()
data = gut.readline().split()
Ngut = int(data[1])
gut.close()

jastrow = open('jastrowidx.def', 'r')
data = jastrow.readline()
data = jastrow.readline().split()
Njastrow = int(data[1])
jastrow.close()

Ndh2, Ndh4, Nidx = 0, 0, NOrbital

fname = 'zqp_opt.dat'
f = open(fname, 'w')
for idx in range(2+Ngut+Njastrow+6*Ndh2+10*Ndh4):
    f.write('{0:25.18e} {1:25.18e} 0.0 '.format(0.0, 0.0))
for idx in range(NOrbital):
    i = orbitalSiteI[idx]
    j = orbitalSiteJ[idx]
    if APFlag == 1:
        isub = subIndex(i)
        dj = positionToIndex(direction(i, j))
        sgn = sgnAP(i, direction(i, j))
        f.write('{0:25.18e} {1:25.18e} 0.0 '.format(fRe[i+Nsite*j]*sgn, 0.0))
    else:
        f.write('{0:25.18e} {1:25.18e} 0.0 '.format(fRe[i+Nsite*j], 0.0))
f.close()
