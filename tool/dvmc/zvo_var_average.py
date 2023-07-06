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

import sys
import numpy as np


def main():
    """Docs
    """
    if (len(sys.argv) == 5):
        fileName1 = sys.argv[1]
        N_ave = float(sys.argv[2])
        average_E = float(sys.argv[3])
        tol_E = float(sys.argv[4])
    else:
        print('example:\n$ zvo_var_average.py zbo_var_001.dat 500'
              ' -101.95 1.0 \n\n')
        sys.exit()

    zvo_var_001 = open(fileName1).read()
    zvo_var_lines = zvo_var_001.split('\n')
    N = len(zvo_var_lines)
    n = len(zvo_var_lines[0].split())
    if n % 3 != 0:
        print("error: incorrect number of elements\n")
        print(n)
        sys.exit()

    zqp_opt = np.zeros(n)
    N_ave_official = 0

    for ii in range(N):
        if ii > N-N_ave:
            if len(zvo_var_lines[ii]) > 0:
                line = zvo_var_lines[ii].split()
                if (abs(average_E - float(line[0])) < tol_E):
                    for jj in range(n):
                        zqp_opt[jj] += float(line[jj])
                    N_ave_official += 1

    zqp_opt = zqp_opt / (float(N_ave_official))
    zqp_opt_file = open("zqp_opt.tmp", 'w')

    for jj in range(n):
        zqp_opt_file.write("%1.18e  " % zqp_opt[jj])
    zqp_opt_file.close()

    return


if __name__ == "__main__":
    main()
