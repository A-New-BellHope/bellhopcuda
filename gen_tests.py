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

dims = {
    2: '2D',
    3: '3D',
    4: 'Nx2D'
}
run_types = {
    'R': 'ray',
    'C': 'coherent',
    'S': 'semicoherent',
    'I': 'incoherent',
    'E': 'eigenrays',
    'A': 'arrivalsascii',
    'a': 'arrivalsbinary'
}
infl_types = {
    'R': 'cervenyraycen',
    'C': 'cervenycart',
    'G': 'hatcartG',
    '^': 'hatcartcaret',
    ' ': 'hatcartspace',
    'g': 'hatraycen',
    'B': 'gaussiancart',
    'b': 'gaussianraycen',
    'S': 'sgb'
}
ssp_types = {
    'N': 'n2linear',
    'C': 'clinear',
    'S': 'cubic',
    'P': 'pchip',
    'Q': 'quad',
    'H': 'hexahedral',
    'A': 'analytic'
}

def should_work(dim, rt, it, st):
    if st in {'P', 'Q'} and dim != 2:
        return False
    if st == 'H' and dim == 2:
        return False
    if it in {'R', 'C', 'S'} and dim == 3:
        return False
    if it == 'C' and dim == 4:
        return False
    if it == 'b' and dim == 2:
        return False
    if it in {'R', 'C'} and rt in {'E', 'A', 'a'}:
        return False
    return True

def gen_field_test(dim, rt, it, st):
    subnames = [dims[dim], run_types[rt], infl_types[it], ssp_types[st]]
    env_name = 'gen_' + '_'.join(subnames)
    print(env_name)
    with open('test/in/' + env_name + '.env', 'w') as envfil:
        envfil.write('\'Gen: ' + ', '.join(subnames) + '\' ! TITLE\n')
        envfil.write('50.0       ! FREQ (Hz)\n')
        envfil.write('1          ! NMEDIA\n')
        envfil.write('\'' + st + 'VW - \'   ! SSP (' + ssp_types[st] 
            + '), top bc (vacuum), atten units, add vol atten, altimetry, dev mode\n')
        NPts = 3
        envfil.write('0  0.0 5000.0 ! NPts (ignored), Sigma (ignored), bot depth\n')
        if st != 'A':
            envfil.write('   0.0 1547.0 /\n')
            envfil.write('1234.5 1500.0 /\n')
            envfil.write('5000.0 1560.0 /\n')
        envfil.write('\'R-    \' 0.0  ! bot bc (rigid), bathymetry, 4 spaces; Sigma (printed but ignored)\n')
        if dim != 2:
            envfil.write('2             ! NSX\n')
            envfil.write('-20.0 20.0  / ! SX(1:NSX) (km)\n')
            envfil.write('2             ! NSY\n')
            envfil.write('-20.0 20.0  / ! SY(1:NSY) (km)\n')
        envfil.write('2             ! NSD\n')
        envfil.write('347.0 682.0 / ! SD(1:NSD) (m)\n')
        if rt in {'E', 'A', 'a'}:
            envfil.write('2               ! NRD\n')
            envfil.write('1135.8 1145.8 / ! RD(1:NRD) (m)\n')
            envfil.write('2               ! NR\n')
            envfil.write('37.2  37.21 /   ! R(1:NR ) (km)\n')
            if dim != 2:
                envfil.write('2             ! Ntheta (number of bearings)\n')
                envfil.write('43.2  43.9 /  ! bearing angles (degrees)\n')
        else:
            envfil.write('51            ! NRD\n')
            envfil.write('0.0 5000.0 /  ! RD(1:NRD) (m)\n')
            envfil.write('51            ! NR\n')
            envfil.write('0.0  100.0 /  ! R(1:NR ) (km)\n')
            if dim != 2:
                envfil.write('10            ! Ntheta (number of bearings)\n')
                envfil.write('0.0  360.0 /  ! bearing angles (degrees)\n')
        envfil.write('\'' + rt + it + ' RR' + ('3' if dim == 3 else '2') 
            + '\'      ! RunType, infl/beam type, ignored, point source, rectilinear grid, dim\n')
        envfil.write('21            ! NBEAMS\n')
        envfil.write('-51.2 51.2 /  ! ALPHA1, 2 (degrees)\n')
        if dim != 2:
            envfil.write('7             ! Nbeta\n')
            envfil.write('0.0 360.0 /   ! beta1, beta2 (degrees) bearing angle fan\n')
        if dim == 2:
            envfil.write('1000.0 5500.0 101.0 ! deltas, box depth, box range\n')
        else:
            envfil.write('1000.0 101.0 101.0 5500.0 ! deltas, box X, box Y, box depth\n')
        if it in {'R', 'C'}:
            envfil.write('\'MS\' 1.0 100.0 0, ! \'Width Curvature\' epsMultiplier rLoop ISINGL (ignored)\n')
            envfil.write('1 4 \'P\' ! Nimage iBeamWindow Component\n')
    if st in {'Q', 'H'}:
        with open('test/in/' + env_name + '.ssp', 'w') as sspfil:
            if st == 'Q':
                Nr, Rmin, Rmax = 5, -102.0, 102.0
                def gen_ssp(r, z):
                    return 1500.0 + 50.0 * r / Nr + 50.0 * z / NPts
                sspfil.write(str(Nr) + '\n')
                sspfil.write(' '.join('{:.2f}'.format(Rmin+(Rmax-Rmin)*r/(Nr-1)) for r in range(Nr)) + '\n')
                for z in range(NPts):
                    sspfil.write(' '.join('{:.2f}'.format(gen_ssp(r, z)) for r in range(Nr)) + '\n')
            else:
                Nx, Ny, Nz, xymin, xymax, zmax = 4, 6, 5, -150.0, 150.0, 5000.0
                def gen_ssp(x, y, z):
                    return 1500.0 + 50.0 * x / Nx - 50.0 * y / Ny + 50.0 * z / Nz
                sspfil.write(str(Nx) + '\n')
                sspfil.write(' '.join('{:.2f}'.format(xymin+(xymax-xymin)*x/(Nx-1)) for x in range(Nx)) + '\n')
                sspfil.write(str(Ny) + '\n')
                sspfil.write(' '.join('{:.2f}'.format(xymin+(xymax-xymin)*y/(Ny-1)) for y in range(Ny)) + '\n')
                sspfil.write(str(Nz) + '\n')
                sspfil.write(' '.join('{:.2f}'.format(zmax*z/(Nz-1)) for z in range(Nz)) + '\n')
                for z in range(Nz):
                    for y in range(Ny):
                        sspfil.write(' '.join('{:.2f}'.format(gen_ssp(x, y, z)) for x in range(Nx)) + '\n')
    return env_name
    
def gen_all_it_st_combos(passtxt, failtxt, dim, rt):
    it_list = ['G'] if rt == 'R' else infl_types.keys()
    for it in it_list:
        for st in ssp_types.keys():
            env_name = gen_field_test(dim, rt, it, st)
            (passtxt if should_work(dim, rt, it, st) else failtxt).write(env_name + '\n')

for dim in {2, 3}:
    dimname = '3d' if dim == 3 else ''
    print('\n\nray{}:\n'.format(dimname))
    with open('gen_ray{}_pass.txt'.format(dimname), 'w') as passtxt, \
        open('gen_ray{}_fail.txt'.format(dimname), 'w') as failtxt:
        gen_all_it_st_combos(passtxt, failtxt, dim, 'R')
    print('\n\ntl{}:\n'.format(dimname))
    with open('gen_tl{}_pass.txt'.format(dimname), 'w') as passtxt, \
        open('gen_tl{}_fail.txt'.format(dimname), 'w') as failtxt:
        for rt in ['C', 'S', 'I']:
            gen_all_it_st_combos(passtxt, failtxt, dim, rt)
    print('\n\neigen{}:\n'.format(dimname))
    with open('gen_eigen{}_pass.txt'.format(dimname), 'w') as passtxt, \
        open('gen_eigen{}_fail.txt'.format(dimname), 'w') as failtxt:
        gen_all_it_st_combos(passtxt, failtxt, dim, 'E')
    print('\n\narr{}:\n'.format(dimname))
    with open('gen_arr{}_pass.txt'.format(dimname), 'w') as passtxt, \
        open('gen_arr{}_fail.txt'.format(dimname), 'w') as failtxt:
        for rt in ['A', 'a']:
            gen_all_it_st_combos(passtxt, failtxt, dim, rt)
