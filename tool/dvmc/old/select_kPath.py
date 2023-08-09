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
import sys

StdFace = open('StdFace.def').read().replace(' ', '')
U = 0.
L = W = 1
if (StdFace.find('L=') >= 0):
    L = int(re.compile('L=([0-9]*)').findall(StdFace)[0])
if (StdFace.find('W=') >= 0):
    W = int(re.compile('W=([0-9]*)').findall(StdFace)[0])
if (StdFace.find('U=') >= 0):
    U = float(re.compile('U=([0-9.]*)').findall(StdFace)[0])


outputDir = 'output/'
spectrumparaFileName = 'spectrumpara.def'

if (len(sys.argv) == 1):
    pass
elif (len(sys.argv) == 2):
    spectrumparaFileName = sys.argv[1]
elif (len(sys.argv) == 3):
    outputDir = sys.argv[2]+'/'
else:
    print('example:\n$ select_kPath.py \nor:\n$ select_kPath.py'
          '  spectrumpara.def\nor:\n$ select_kPath.py spectrumpara.def output/')
    sys.exit()


def main():

    spectrumpara = open(spectrumparaFileName).read()

    for line in spectrumpara.split('\n'):
        if len(line) > 0:
            if line[0] != '#':
                term = line.split()
                if len(term) > 0:
                    if (term[0][:] == 'kPath' or term[0][:-1] == 'kPath'):
                        if term[1] == 'all':
                            kPath = range(W*L)
                        elif term[1][0:6] == 'range(':
                            kPath = range(int(term[1][6:-1]))
                        else:
                            kPath = ReadRange(term[1])

    def select_kPath(fileName, fileOut):
        Akw_all = open(fileName).read()
        Akw = open(fileOut, 'w')
        for line in Akw_all.split('\n'):
            if len(line) > 0:
                terms = line.split()
                lineToPrint = ''
                for k in kPath:
                    lineToPrint += terms[int(k)] + '  '
            Akw.write(lineToPrint+'\n')
        Akw.close()

    select_kPath(outputDir+'Akw_all.dat',   outputDir+'Akw.dat')
    select_kPath(outputDir+'Akw_e_all.dat', outputDir+'Akw_e.dat')
    select_kPath(outputDir+'Akw_h_all.dat', outputDir+'Akw_h.dat')

    return


def ReadRange(inputStr: str) -> list:
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
            print('sepated by comma. Pairs must contains only one ":"')
            print('example:')
            print('0:3,5,6,8:11')
            print('will be interpreted as')
            print('[0,1,2,3,5,6,8,9,10,11]')
            exit()

    return rangeOut


if __name__ == "__main__":
    main()
