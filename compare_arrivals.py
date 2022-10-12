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

def compare_asc_files(cxxf, forf):
    l = 0
    cxxmorelines = 0
    errored = False
    while True:
        l += 1
        cxxl = cxxf.readline()
        forl = forf.readline()
        if not cxxl and not forl:
            break
        def fatalerror():
            print('Parse error: line {}:'.format(l))
            print('FOR: ' + forls)
            print('CXX: ' + cxxls)
            sys.exit(1)
        def error():
            print('Results failed to match: line {}:'.format(l))
            print('FOR: ' + forls)
            print('CXX: ' + cxxls + '\n')
            nonlocal errored
            errored = True
        def findtype(t):
            return 'str' if '\'' in t else 'float' if '.' in t else 'int'
        cxxtokens = cxxl.split()
        fortokens = forl.split()
        cxxls = ' '.join(cxxtokens)
        forls = ' '.join(fortokens)
        if len(cxxtokens) != len(fortokens): fatalerror()
        if len(cxxtokens) == 1 and l >= 7:
            if cxxmorelines != 0:
                print('Number of extra lines didn\'t match actual number')
                sys.exit(1)
            cxxmorelines = int(cxxtokens[0]) - int(fortokens[0])
            if cxxmorelines != 0:
                print('Line {}: FOR {} arrivals / CXX {} arrivals'.format(l, fortokens[0], cxxtokens[0]))
        else:
            t = 0
            for cxxt, fort in zip(cxxtokens, fortokens):
                ty = findtype(fort)
                if ty != findtype(cxxt): fatalerror()
                if ty == 'str':
                    if cxxt != fort: fatalerror()
                elif ty == 'int':
                    if int(cxxt) != int(fort): fatalerror()
                else:
                    cf = float(cxxt)
                    ff = float(fort)
                    if ff == 0.0:
                        if cf != 0.0:
                            error()
                            break
                        continue
                    if abs((cf - ff) / ff) <= 1e-5: continue
                    if t != 0:
                        error()
                        break
                    if cxxmorelines > 0:
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
                        cxxmorelines -= 1
                    elif cxxmorelines < 0:
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
                        cxxmorelines += 1
                        l += 1
                    else:
                        error()
                    break
                t += 1
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
