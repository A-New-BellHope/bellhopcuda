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

use_float = False
if use_float:
    thresh_rel_print = 1e-4
    thresh_rel_alarm = 1e-2
else:
    thresh_rel_print = 1e-6
    thresh_rel_alarm = 1e-3

def compare_files(cxxf, forf):
    def compare_lines(cxxl, forl):
        def printlns():
            print('CXX: {}\nFOR: {}'.format(cxxl, forl))
        if "'" in cxxl:
            if not "'" in forl:
                return False
            if cxxl != forl:
                printlns()
            return True
        cxxt = [t for t in cxxl.split() if t]
        fort = [t for t in forl.split() if t]
        if len(cxxt) != len(fort):
            return False
        for i, (c, f) in enumerate(zip(cxxt, fort)):
            if '.' in c:
                if '.' not in f:
                    return False
                try:
                    cf = float(c)
                    ff = float(f)
                except ValueError:
                    return False
                if abs(cf - ff) < 1e-8:
                    continue # Ignore relative error if abs error very small
                if (cf == 0.0) != (ff == 0.0):
                    return False
                if cf != 0.0:
                    relerror = abs((cf - ff) / ff)
                    if relerror > thresh_rel_alarm:
                        return False
                    elif relerror > thresh_rel_print:
                        printlns()
            else:
                try:
                    ci = int(c)
                    fi = int(f)
                except ValueError:
                    return False
                if ci == fi:
                    continue
                elif len(cxxt) == 3 and i == 0 and \
                    all('.' not in t for t in (cxxt + fort)) and \
                    int(cxxt[1]) == int(fort[1]) and int(cxxt[2]) == int(fort[2]):
                    printlns()
                    if abs(ci - fi) > 1:
                        print('Warning, CXX {} steps / FOR {} steps'.format(ci, fi))
                else:
                    return False
        return True
    cxx_read = True
    for_read = True
    cxx_num = 0
    for_num = 0
    while True:
        if cxx_read:
            cxxl = cxxf.readline().strip()
            cxx_num += 1
        if for_read:
            forl = forf.readline().strip()
            for_num += 1
        if not cxxl and not forl:
            break
        cxx_read = for_read = True
        res = compare_lines(cxxl, forl)
        if not res:
            cxxl2 = cxxf.readline().strip()
            forl2 = forf.readline().strip()
            if compare_lines(cxxl, forl2):
                print('Extra FOR line {}: {}'.format(for_num, forl))
                cxx_read = False
                cxxl = cxxl2
            elif compare_lines(cxxl2, forl):
                print('Extra CXX line {}: {}'.format(cxx_num, cxxl))
                for_read = False
                forl = forl2
            else:
                print('Ray files failed to match\nCXX line {}: {}\nFOR line {}: {}'.format(
                    cxx_num, cxxl, for_num, forl))
                sys.exit(1)
            cxx_num += 1
            for_num += 1

if len(sys.argv) not in {2, 3}:
    print('Usage: python3 compare_ray.py MunkB_ray [cxx1/cxxmulti/cuda]')
    print('No paths, no .ray')
    sys.exit(1)

rayfil = sys.argv[1]
if len(sys.argv) == 3:
    comparisons = [sys.argv[2]]
else:
    comparisons = ['cxx1', 'cxxmulti', 'cuda']

for c in comparisons:
    with open('test/FORTRAN/{}.ray'.format(rayfil), 'r') as forf:
        cxxfile = 'test/{}/{}.ray'.format(c, rayfil)
        try:
            with open(cxxfile, 'r') as cxxf:
                print('Ray comparison FORTRAN vs. {}:'.format(c))
                compare_files(cxxf, forf)
        except FileNotFoundError:
            print('{} not found, skipping {}'.format(cxxfile, c))
