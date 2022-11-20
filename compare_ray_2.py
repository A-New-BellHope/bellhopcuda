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

print_matches = False
use_float = False
if use_float:
    thresh_rel_print = 1e-4
    thresh_rel_alarm = 1e-2
else:
    thresh_rel_print = 1e-6
    thresh_rel_alarm = 1e-3

def load_rayfil(f):
    ret = {'rays': []}
    ray = None
    state = 0
    n2 = None
    for i, l in enumerate(f):
        l = l.strip()
        t = [x for x in l.split() if x]
        if state == 0:
            assert l.startswith('\'') and l.endswith('\'')
            ret['Title'] = l[1:-1]
            state = 1
        elif state == 1:
            ret['freq0'] = float(l)
            state = 2
        elif state == 2:
            assert len(t) == 3
            ret['NSx'], ret['NSy'], ret['NSz'] = int(t[0]), int(t[1]), int(t[2])
            state = 3
        elif state == 3:
            assert len(t) == 2
            ret['Nalpha'], ret['Nbeta'] = int(t[0]), int(t[1])
            state = 4
        elif state == 4:
            ret['Top.Depth'] = float(l)
            state = 5
        elif state == 5:
            ret['Bot.Depth'] = float(l)
            state = 6
        elif state == 6:
            if l == '\'rz\'':
                ret['3D'] = False
            elif l == '\'xyz\'':
                ret['3D'] = True
            else:
                raise RuntimeError
            state = 7
        else:
            if state == 9 and n2 == 0:
                ret['rays'].append(ray)
                state = 7
            if state == 7:
                ray = {'alpha0': float(l), 'line': i+1, 'data': []}
                state = 8
            elif state == 8:
                assert len(t) == 3
                n2, ray['NumTopBnc'], ray['NumBotBnc'] = int(t[0]), int(t[1]), int(t[2])
                state = 9
            elif state == 9:
                assert len(t) == 2
                ray['data'].append((float(t[0]), float(t[1])))
                n2 -= 1
            else:
                raise RuntimeError
    assert state == 9 and n2 == 0
    ret['rays'].append(ray)
    return ret

def compare_floats(cf, ff, cl, fl):
    if abs(cf - ff) < 1e-8:
        return True # Ignore relative error if abs error very small
    if (cf == 0.0) != (ff == 0.0):
        return False
    if cf != 0.0:
        relerror = abs((cf - ff) / ff)
        if relerror > thresh_rel_alarm:
            return False
        elif relerror > thresh_rel_print:
            print('cxx line {:7} / for line {:7}: for value {:24.17f} rel err {}'.format(
                cl, fl, ff, relerror))
    return True

def compare_points(cp, fp, cl, fl):
    return all(compare_floats(cp[i], fp[i], cl, fl) for i in range(2))

def compare_rays(cxxr, forr, exact):
    if any(cxxr[k] != forr[k] for k in {'NumTopBnc', 'NumBotBnc'}) or not compare_floats(
        cxxr['alpha0'], forr['alpha0'], cxxr['line'], forr['line']):
        # print('cxx and for rays had different metadata, lines {} and {}'.format(
        #     cxxr['line'], forr['line']))
        return False
    cd, fd = cxxr['data'], forr['data']
    ci, fi = 0, 0
    while True:
        cl, fl = cxxr['line'] + ci, forr['line'] + fi
        if (ci >= len(cd)) != (fi >= len(fd)):
            print('cxx and for rays ended at different steps, lines {} and {}'.format(
                cl, fl))
            return False
        if ci >= len(cd):
            return True
        if not compare_points(cd[ci], fd[fi], cl, fl):
            if not exact and ci < len(cd) - 1 and compare_points(cd[ci+1], fd[fi], cl+1, fl):
                print('Extra cxx step line {}: {}'.format(cl, cd[ci]))
                ci += 1
            elif not exact and fi < len(fd) - 1 and compare_points(cd[ci], fd[fi+1], cl, fl+1):
                print('Extra for step line {}: {}'.format(fl, fd[fi]))
                fi += 1
            else:
                if not exact:
                    print('cxx line {} did not match for line {}:\ncxx {}\nfor {}'.format(
                        cl, fl, cd[ci], fd[fi]))
                return False
        ci += 1
        fi += 1

def find_matching_ray(cxxr, ford, exact):
    fori = -1
    while True:
        fori += 1
        if fori >= len(ford['rays']):
            break
        forr = ford['rays'][fori]
        if compare_rays(cxxr, forr, exact):
            break
    if fori < len(ford['rays']):
        if print_matches:
            print('cxx ray line {} matched for ray line {}'.format(cxxr['line'], forr['line']))
        del ford['rays'][fori]
        return True
    return False

def compare_raydata(cxxd, ford):
    for k in {'freq0', 'NSx', 'NSy', 'NSz', 'Nalpha', 'Nbeta', 'Top.Depth',
        'Bot.Depth', '3D'}:
        if cxxd[k] != ford[k]:
            print('Ray files failed to match: {} cxx {} for {}'.format(k, cxxd[k], ford[k]))
            sys.exit(1)
    if len(cxxd['rays']) != len(ford['rays']):
        print('Ray files failed to match: cxx {} rays, for {} rays'.format(
            len(cxxd['rays']), len(ford['rays'])))
        sys.exit(1)
    for cxxr in cxxd['rays']:
        if find_matching_ray(cxxr, ford, True):
            continue
        if find_matching_ray(cxxr, ford, False):
            continue
        print('Ray files failed to match: no for ray found which matches '
            + 'cxx ray starting on line {}'.format(cxxr['line']))
        sys.exit(1)

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
        ford = load_rayfil(forf)
    cxxfile = 'test/{}/{}.ray'.format(c, rayfil)
    try:
        with open(cxxfile, 'r') as cxxf:
            print('Ray comparison FORTRAN vs. {}:'.format(c))
            cxxd = load_rayfil(cxxf)
        compare_raydata(cxxd, ford)
    except FileNotFoundError:
        print('{} not found, skipping {}'.format(cxxfile, c))
