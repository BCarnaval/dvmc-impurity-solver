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

if len(sys.argv) != 2:
    print("example:\n$ python cluster.py 12")
    sys.exit()

NN = int(sys.argv[1])

for ii in range(NN):
    print('#')
    for jj in range(NN):
        print(' %3d' % (NN * NN - (ii + 1) * NN + jj))
    print('')
