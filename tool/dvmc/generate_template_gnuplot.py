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

full_path = os.path.realpath(__file__)
pythonPathCode, file1 = os.path.split(full_path)

if os.path.isfile('StdFace.def'):
    StdFace = open('StdFace.def').read().replace(' ', '')
    if (StdFace.find('L=') >= 0):
        L = int(re.compile('L=([0-9]*)').findall(StdFace)[0])
    if (StdFace.find('W=') >= 0):
        W = int(re.compile('W=([0-9]*)').findall(StdFace)[0])
    if (StdFace.find('U=') >= 0):
        U = float(re.compile('U=([0-9.]*)').findall(StdFace)[0])
elif os.path.isfile('params'):
    paramfile = open('params', 'r')
    for line in paramfile:
        param_name = line.split()[0]
        param_val = line.split()[1]
        if (param_name == 'L'):
            L = int(param_val)
        if (param_name == 'W'):
            W = int(param_val)
        if (param_name == 'U'):
            U = float(param_val)
else:
    print("no input file found!")
    sys.exit()


spectrumparaFileName = 'spectrumpara.def'

if (len(sys.argv) == 1):
    print('')
elif (len(sys.argv) == 2):
    spectrumparaFileName = sys.argv[1]
else:
    print("example:\n$ vmc_spectrum.py \nor:\n$ vmc_spectrum.py spectrumpara.def")
    sys.exit()


def main() -> None:
    """Docs
    """
    w_min_data, w_max_data = -15.0, 15.0
    w_min, w_max = -8.0, 8.0
    spectrumpara = open(spectrumparaFileName).read()

    for line in spectrumpara.split('\n'):
        if len(line) > 0:
            if line[0] != '#':
                term = line.split()
                if len(term) > 0:
                    if term[0] == 'w_min_data':
                        w_min_data = float(term[1])

                    if term[0] == 'w_max_data':
                        w_max_data = float(term[1])

                    if term[0] == 'w_min':
                        w_min = float(term[1])

                    if term[0] == 'w_max':
                        w_max = float(term[1])

                    if term[0] == 'Nw':
                        Nw = int(term[1])

                    if term[0][:] == 'kPath' or term[0][:-1] == 'kPath':
                        if term[1] == 'all':
                            kPath = range(W * L)

                        elif term[1][0:6] == 'range(':
                            kPath = range(int(term[1][6:-1]))

                        else:
                            kPath = ReadRange(term[1])

                        if term[0][-1] == '1':
                            kPath1 = kPath

                        elif term[0][-1] == '2':
                            pass

                        else:
                            kPath1 = kPath

    def replaceTemplateFileValues(templateFileName: str) -> None:
        """Docs
        """
        gnuplotString = open(pythonPathCode + '/' + templateFileName).read()
        gnuplotString = gnuplotString.replace(
            'w_min_data = -13.0', 'w_min_data =% 4.5f' % (w_min_data))
        gnuplotString = gnuplotString.replace(
            'w_max_data =  13.0', 'w_max_data =% 4.5f' % (w_max_data))
        gnuplotString = gnuplotString.replace(
            'w_min = -8.0', 'w_min =% 4.5f' % (w_min))
        gnuplotString = gnuplotString.replace(
            'w_max =  8.0', 'w_max =% 4.5f' % (w_max))
        gnuplotString = gnuplotString.replace(
            'kRange = 32', 'kRange =% 4.5f' % (len(kPath1)-1))
        gnuplotString = gnuplotString.replace('Nw = 1500', 'Nw =% 4.5f' % (Nw))

        fileOut = open('./' + templateFileName[9:], 'w')
        fileOut.write(gnuplotString)
        fileOut.close()

        return

    replaceTemplateFileValues('template_plot_allAkw.gp')
    replaceTemplateFileValues('template_plot_Akw.gp')
    replaceTemplateFileValues('template_plot_dos.gp')

    return

###############################################################################
# Sub routines
###############################################################################


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
            print('sepated by comma. Pairs must contains only one ":"')
            print('example:')
            print('0:3,5,6,8:11')
            print('will be interpreted as')
            print('[0,1,2,3,5,6,8,9,10,11]')
            exit()

    return rangeOut


if __name__ == "__main__":
    main()
