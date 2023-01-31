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

def parsems(s):
    a = s.split('m')
    assert len(a) == 2
    assert a[1].endswith('s')
    m = int(a[0])
    s = float(a[1][:-1])
    t = (60 * m + s) * 1000.0
    return str(round(t))
    
def parsecxx(s):
    f = float(s)
    if f >= 1000.0:
        return str(int(round(f, 0)))
    else:
        return str(round(f, 1))

with open(sys.argv[1], 'r') as o:
    f = None
    runs = []
    writes = []
    def printvals():
        if f is not None:
            assert len(runs) == len(writes) == 3
            print('\t'.join([f] + runs + writes))
    for l in o:
        l = l.strip()
        toks = l.split()
        if l.startswith('genperf_'):
            printvals()
            f = None
            runs = []
            writes = []
        elif l.startswith('real\t'):
            f = parsems(toks[1])
        elif l.startswith('Run: '):
            runs.append(parsecxx(toks[1]))
        elif l.startswith('FinalizeXMode: '):
            writes.append(parsecxx(toks[1]))
    printvals()
