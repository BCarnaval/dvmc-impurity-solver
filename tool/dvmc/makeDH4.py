#!/usr/bin/env python3
import sys
import re

StdFace = open('StdFace.def').read().replace(' ', '')
L = int(re.compile('L=([0-9]*)').findall(StdFace)[0])
W = int(re.compile('W=([0-9]*)').findall(StdFace)[0])


def position(r):
    """Docs
    """
    x = (r[0]+W) % W
    y = (r[1]+L) % L
    return (x, y)


def indexToPosition(i):
    """Docs
    """
    x = i % W
    y = i/W
    return (x, y)


def positionToIndex(r):
    """Docs
    """
    x = (r[0]+W) % W
    y = (r[1]+L) % L
    return x+W*y


def neighborIndex(i, dr):
    """Docs
    """
    r = indexToPosition(i)
    x = r[0]+dr[0]
    y = r[1]+dr[1]
    return positionToIndex([x, y])


if len(sys.argv) > 1:
    print("Make_DH.py")
    sys.exit()

Nsite = L*W


# DH4.def
NDoublonHolon4 = 1
f = open('DH4.def', 'w')
f.write('==============================\n' +
        "NDoublonHolon4SiteIdx\t" + str(NDoublonHolon4) + "\n" +
        "ComplexType  0\n" +
        '==============================\n' +
        '==============================\n')

t = 0  # nn_+
for i in range(Nsite):
    j1 = neighborIndex(i, [1, 0])
    j2 = neighborIndex(i, [-1, 0])
    j3 = neighborIndex(i, [0, 1])
    j4 = neighborIndex(i, [0, -1])
    f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(i, j1, j2, j3, j4, t))

activate = 1
for i in range(10*NDoublonHolon4):
    f.write("{0}\t{1}\n".format(i, activate))
f.close()
