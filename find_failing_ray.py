'''
bellhopcxx / bellhopcuda - C++/CUDA port of BELLHOP underwater acoustics simulator
Copyright (C) 2021-2022 The Regents of the University of California
c/o Jules Jaffe team at SIO / UCSD, jjaffe@ucsd.edu
Based on BELLHOP, which is Copyright (C) 1983-2020 Michael B. Porter

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see <https://www.gnu.org/licenses/>.
'''
import sys

if sys.version_info.major < 3:
    print('This is a python3 script')
    sys.exit(-1)

import random
import os
import shutil

if len(sys.argv) != 3:
    print('Usage: python3 find_failing_ray.py (ray/tl)(2D/3D/Nx2D) munk_something')
    print('No path, no .env')
    sys.exit(1)

testpath = 'test/in/'
envfil = sys.argv[2]
envfilename = testpath + envfil + '.env'
if 'Nx2D' in sys.argv[1]:
    dim = 4
elif '3D' in sys.argv[1]:
    dim = 3
elif '2D' in sys.argv[1]:
    dim = 2
else:
    raise RuntimeError('Invalid run type')
if 'ray' in sys.argv[1]:
    compare_script = 'compare_ray_2.py'
elif 'tl' in sys.argv[1]:
    compare_script = 'compare_shdfil.py'
else:
    raise RuntimeError('Invalid run type')
alphathreesixtysweep = False
betathreesixtysweep = True

def set_random_ray():
    tempfilename = 'temp12345.env'
    alpha = beta = None
    with open(envfilename, 'r') as infile, open(tempfilename, 'w') as tempfile:
        for l in infile:
            isalpha = any(k in l.lower() for k in {'nalpha', 'nbeams'})
            isbeta = 'nbeta' in l.lower()
            assert not (isalpha and isbeta)
            if isalpha or isbeta:
                if isalpha and alpha is not None:
                    print('Already modified alpha, got "' + l + '"')
                    sys.exit(2)
                if isbeta and beta is not None:
                    print('Already modified beta, got "' + l + '"')
                    sys.exit(2)
                halves = l.split('!')
                assert len(halves) == 2
                toks = halves[0].split()
                if len(toks) == 1:
                    print('Need Nalpha and iSingleAlpha (or same for beta) on one line, got "' + l + '"')
                    sys.exit(2)
                assert len(toks) == 2
                n = int(toks[0])
                assert n > 1
                threesixtysweep = alphathreesixtysweep if isalpha else betathreesixtysweep
                i = random.randint(1, n-1 if isbeta and threesixtysweep else n)
                if isalpha:
                    alpha = i
                else:
                    beta = i
                l = str(n) + ' ' + str(i) + '             !' + halves[1]
            tempfile.write(l)
    if alpha is None:
        print('Did not find Nalpha line')
        sys.exit(2)
    os.replace(tempfilename, envfilename)
    return alpha, beta

while True:
    alpha, beta = set_random_ray()
    shutil.copyfile(envfilename, 'test/cxx1/' + envfil + '.env')
    shutil.copyfile(envfilename, 'test/FORTRAN/' + envfil + '.env')
    if os.system('./bin/bellhopcxx -1 -' + str(dim) + ' test/cxx1/' + envfil) != 0:
        print('bellhopcxx failed')
        sys.exit(3)
    if os.system('../bellhop/Bellhop/bellhop' + ('' if dim == 2 else '3d') + '.exe test/FORTRAN/' + envfil) != 0:
        print('BELLHOP failed')
        sys.exit(4)
    passed = os.system('python3 ' + compare_script + ' ' + envfil + ' cxx1') == 0
    print('{} with alpha {} beta {}'.format('passed' if passed else 'failed', alpha, beta))
    if not passed:
        break
