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

import struct
import difflib
from math import isfinite

def compare_asc_files(cxxf, forf):
    def fatalerror(msg = 'Parse error'):
        print('{}: line {}:'.format(msg, l))
        print('FOR: ' + forls)
        print('CXX: ' + cxxls)
        raise RuntimeError
    def error():
        print('Results failed to match: line {}:'.format(l))
        print('FOR: ' + forls)
        print('CXX: ' + cxxls + '\n')
        nonlocal errored
        errored = True
    def findtype(t):
        if '\'' in t:
            return 'str'
        if any(x in t for x in {'.', 'nan', 'inf'}):
            return 'float'
        if all(x.isnumeric() or x == '-' for x in t):
            return 'int'
        fatalerror()
    def compare_floats(cxxt, fort):
        cf = float(cxxt)
        ff = float(fort)
        if not isfinite(cf) or not isfinite(ff):
            fatalerror('Non-finite results')
        if abs(cf - ff) < 1e-8: return True
        if (ff == 0.0) != (cf == 0.0): return False
        if ff == 0.0: return True
        if abs((cf - ff) / ff) <= 1e-5: return True
        return False
    def compare_tokens(cxxtokens, fortokens):
        for cxxt, fort in zip(cxxtokens, fortokens):
            ty = findtype(fort)
            if ty != findtype(cxxt): return False
            if ty == 'str':
                if cxxt != fort: return False
            elif ty == 'int':
                if int(cxxt) != int(fort): return False
            else:
                return compare_floats(cxxt, fort)
        return True
    l = 0
    cxxmorelines = 0
    errored = False
    threed = None
    datamode = False
    maxline = False
    countline = False
    endoffile = False
    NSx, NSy, NSz, NRz, NRr, Ntheta = 1, 1, None, None, None, 1
    sx, sy, sz, rz, rr, theta = 0, 0, 0, 0, 0, 0
    while True:
        l += 1
        cxxl = cxxf.readline()
        forl = forf.readline()
        if not cxxl and not forl:
            break
        if endoffile:
            fatalerror('Extra lines at end of file')
        cxxtokens = cxxl.split()
        fortokens = forl.split()
        cxxls = ' '.join(cxxtokens)
        forls = ' '.join(fortokens)
        if len(cxxtokens) != len(fortokens): fatalerror()
        nextrcvr = False
        if l == 1:
            if len(cxxtokens) != 1 or cxxtokens[0] != fortokens[0] or \
                cxxtokens[0] not in {'\'2D\'', '\'3D\''}:
                fatalerror()
            threed = cxxtokens[0] == '\'3D\''
        elif not datamode:
            if not compare_tokens(cxxtokens, fortokens):
                fatalerror()
            if l == 2:
                assert len(cxxtokens) == 1 and findtype(cxxtokens[0]) == 'float'
            else:
                assert len(cxxtokens) >= 2 and findtype(cxxtokens[0]) == 'int'
                assert all(findtype(t) == 'float' for t in cxxtokens[1:])
                i = int(cxxtokens[0])
                assert i + 1 == len(cxxtokens)
                if l == 3:
                    if threed:
                        NSx = i
                    else:
                        NSz = i
                elif l == 4:
                    if threed:
                        NSy = i
                    else:
                        NRz = i
                elif l == 5:
                    if threed:
                        NSz = i
                    else:
                        NRr = i
                elif l == 6:
                    assert threed
                    NRz = i
                elif l == 7:
                    assert threed
                    NRr = i
                elif l == 8:
                    assert threed
                    Ntheta = i
            if threed and l == 8 or not threed and l == 5:
                datamode = True
                maxline = True
        elif maxline:
            if len(cxxtokens) != 1 or findtype(cxxtokens[0]) != 'int' or findtype(fortokens[0]) != 'int':
                fatalerror()
            if int(cxxtokens[0]) != int(fortokens[0]):
                print('Max arrivals: FOR {} CXX {}'.format(int(fortokens[0]), int(cxxtokens[0])))
            maxline = False
            countline = True
        elif countline:
            if len(cxxtokens) != 1 or findtype(cxxtokens[0]) != 'int' or findtype(fortokens[0]) != 'int':
                fatalerror()
            cxxarrcount, forarrcount = int(cxxtokens[0]), int(fortokens[0])
            if (cxxarrcount == 0) != (forarrcount == 0):
                fatalerror('One file has no arrivals at this rcvr')
            if cxxarrcount == 0:
                nextrcvr = True
            else:
                if cxxarrcount != forarrcount:
                    print('Line {}: FOR {} arrivals / CXX {} arrivals'.format(l, forarrcount, cxxarrcount))
                countline = False
        else:
            ntok = 10 if threed else 8
            if len(cxxtokens) != ntok: fatalerror()
            for t in range(ntok):
                cxxt = cxxtokens[t]
                fort = fortokens[t]
                ty = findtype(fort)
                if ty != findtype(cxxt): fatalerror()
                if t >= ntok - 2:
                    assert ty == 'int'
                    if int(cxxt) != int(fort):
                        error()
                else:
                    assert ty == 'float'
                    if compare_floats(cxxt, fort):
                        continue
                    if t != 0:
                        error()
                        break
                    cf = float(cxxt)
                    ff = float(fort)
                    if cxxarrcount > forarrcount:
                        cxxnextl = cxxf.readline()
                        cxxtokens = cxxnextl.split()
                        cxxnextls = ' '.join(cxxtokens)
                        if len(cxxtokens) != len(fortokens): fatalerror()
                        cf2 = float(cxxtokens[0])
                        if abs((cf2 + cf - ff) / ff) > 1e-5: fatalerror()
                        print('FOR line corresponds to the sum of these two CXX lines')
                        print('FOR line ' + forls)
                        print('CXX line ' + cxxls)
                        print('CXX line ' + cxxnextls + '\n')
                        cxxarrcount -= 1
                    elif cxxarrcount < forarrcount:
                        fornextl = forf.readline()
                        fortokens = fornextl.split()
                        fornextls = ' '.join(fortokens)
                        if len(fortokens) != len(cxxtokens): fatalerror()
                        ff2 = float(fortokens[0])
                        if abs((cf - (ff + ff2)) / (ff + ff2)) > 1e-5: fatalerror()
                        print('CXX line corresponds to the sum of these two FOR lines')
                        print('CXX line ' + cxxls)
                        print('FOR line ' + forls)
                        print('FOR line ' + fornextls + '\n')
                        forarrcount -= 1
                        l += 1
                    else:
                        error()
                        break
            cxxarrcount -= 1
            forarrcount -= 1
            nextrcvr = True
        if nextrcvr:
            if (cxxarrcount == 0) != (forarrcount == 0):
                fatalerror()
            if cxxarrcount == 0:
                countline = True
                rr += 1
                if rr == NRr:
                    rr = 0
                    rz += 1
                    if rz == NRz:
                        rz = 0
                        theta += 1
                        if theta == Ntheta:
                            theta = 0
                            maxline = True
                            sy += 1
                            if sy == NSy:
                                sy = 0
                                sx += 1
                                if sx == NSx:
                                    sx = 0
                                    sz += 1
                                    if sz == NSz:
                                        endoffile = True
    if not endoffile:
        if not threed and NSz == 1 and NRr == NRz and rr == 0 and rz == 1:
            print('Assuming this is an irregular grid')
        else:
            print(f'Ran out of lines at end of file, source {sx} {sy} {sz} rcvr th_r_z {theta} {rr} {rz}')
            errored = True
    if errored:
        print('Error(s) detected in ASCII arrivals results')
        sys.exit(1)
    else:
        print('ASCII arrival results matched')

def compare_bin_files(cxxf, forf):
    cxxd, ford = cxxf.read(), forf.read()
    if len(cxxd) != len(ford):
        print('CXX and FOR binary files of different length, comparing them is not currently implemented')
        sys.exit(1)
    for a in range(0, len(cxxd), 4):
        ci = struct.unpack('=I', cxxd[a:a+4])[0]
        fi = struct.unpack('=I', ford[a:a+4])[0]
        cf = None
        if fi >= 0x26000000 and fi < 0x5B000000 or fi >= 0xA6000000 and fi < 0xDB000000:
            # Probably a float
            cf = struct.unpack('=f', cxxd[a:a+4])[0]
            ff = struct.unpack('=f', ford[a:a+4])[0]
            if abs((cf - ff) / ff) <= 1e-5: continue
        else:
            # Probably an int
            if ci == fi: continue
        print('At {:08X}, CXX {} FOR {}'.format(a, cxxd[a:a+4].hex(), ford[a:a+4].hex()))
        if cf is not None:
            print('({} vs. {})'.format(cf, ff))
        sys.exit(1)
    print('Binary arrival results matched')

if len(sys.argv) not in {2, 3}:
    print('Usage: python3 compare_arrivals.py MunkB_Arr [cxx1/cxxmulti/cuda]')
    print('No paths, no .arr')
    sys.exit(1)

arrfil = sys.argv[1]
if len(sys.argv) == 3:
    comparisons = [sys.argv[2]]
else:
    comparisons = ['cxx1', 'cxxmulti', 'cuda']

with open('test/FORTRAN/{}.arr'.format(arrfil), 'rb') as forf:
    if forf.read(4) == b'\x04\x00\x00\x00':
        compare_func = compare_bin_files
        open_flag = 'rb'
    else:
        compare_func = compare_asc_files
        open_flag = 'r'

for c in comparisons:
    with open('test/FORTRAN/{}.arr'.format(arrfil), open_flag) as forf:
        cxxfile = 'test/{}/{}.arr'.format(c, arrfil)
        try:
            with open(cxxfile, open_flag) as cxxf:
                print('Arrivals comparison FORTRAN vs. {}:'.format(c))
                compare_func(cxxf, forf)
        except FileNotFoundError:
            print('{} not found, skipping {}'.format(cxxfile, c))
