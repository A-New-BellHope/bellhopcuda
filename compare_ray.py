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

if len(sys.argv) != 2:
    print('Usage: python3 compare_ray.py MunkB_ray')
    print('No path, no .ray')
    sys.exit(1)

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
                elif len(cxxt) == 3 and i == 0 and (ci == fi + 1 or ci + 1 == fi):
                    printlns()
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
            cxx_num += 1
            for_num += 1
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

cxx1file = 'test/cxx1/{}.ray'.format(sys.argv[1])
forfile = 'test/FORTRAN/{}.ray'.format(sys.argv[1])

with open(cxx1file, 'r') as cxxf, open(forfile, 'r') as forf:
    print('bellhopcxx single-threaded')
    compare_files(cxxf, forf)

print('Skipping ray results comparison for bellhopcxx multi-threaded')
print('Skipping ray results comparison for bellhopcuda')
