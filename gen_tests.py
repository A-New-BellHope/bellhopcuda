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
dim_ids = {
    2: ' ',
    3: '3',
    4: '2'
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
    if it == 'b' and dim != 3:
        return False
    if it in {'R', 'C'} and rt in {'E', 'A', 'a'}:
        return False
    return True

def write_env_etc(dim, rt, it, st, p, env_name, title):
    print(env_name)
    with open('test/in/' + env_name + '.env', 'w') as envfil:
        envfil.write('\'Gen: ' + title + '\' ! TITLE\n')
        envfil.write('50.0       ! FREQ (Hz)\n')
        envfil.write('1          ! NMEDIA\n')
        envfil.write('\'' + st + 'VW - \'   ! SSP (' + ssp_types[st] 
            + '), top bc (vacuum), atten units, add vol atten, altimetry, dev mode\n')
        envfil.write('0  0.0 5000.0 ! NPts (ignored), Sigma (ignored), bot depth\n')
        if st != 'A':
            if p['ssp']['NPts'] == 3:
                envfil.write('   0.0 1547.0 /\n')
                envfil.write('1234.5 1500.0 /\n')
                envfil.write('5000.0 1560.0 /\n')
            else:
                assert p['ssp']['NPts'] > 3
                for d in range(p['ssp']['NPts']):
                    z = 5000.0 * d / (p['ssp']['NPts'] - 1)
                    c = max(
                        1547.0 - 47.0 * z / 1234.5,
                        1500.0 + (z - 1234.5) / (5000.0 - 1234.5) * 60.0)
                    envfil.write('{:6.1f} {:6.1f} /\n'.format(z, c))
        envfil.write('\'R-    \' 0.0  ! bot bc (rigid), bathymetry, 4 spaces; Sigma (printed but ignored)\n')
        if dim != 2:
            envfil.write(f"{p['NSx']}            ! NSX\n")
            envfil.write(f"{' '.join(map(str, p['Sx']))} /  ! SX(1:NSX) (km)\n")
            envfil.write(f"{p['NSy']}            ! NSY\n")
            envfil.write(f"{' '.join(map(str, p['Sy']))} /  ! SY(1:NSY) (km)\n")
        envfil.write(f"{p['NSz']}                ! NSD\n")
        envfil.write(f"{' '.join(map(str, p['Sz']))}     /  ! SD(1:NSD) (m)\n")
        envfil.write(f"{p['NRz']}                ! NRD\n")
        envfil.write(f"{' '.join(map(str, p['Rz']))}     /  ! RD(1:NRD) (m)\n")
        envfil.write(f"{p['NRr']}                ! NR\n")
        envfil.write(f"{' '.join(map(str, p['Rr']))}     /  ! R(1:NR ) (km)\n")
        if dim != 2:
            envfil.write(f"{p['Ntheta']}         ! Ntheta (number of bearings)\n")
            envfil.write(f"{' '.join(map(str, p['theta']))} /  ! bearing angles (degrees)\n")
        envfil.write('\'' + rt + it + ' RR' + dim_ids[dim]
            + '\'      ! RunType, infl/beam type, ignored, point source, rectilinear grid, dim\n')
        envfil.write(f"{p['Nalpha']}            ! NBEAMS\n")
        envfil.write(f"{' '.join(map(str, p['alpha']))} /  ! ALPHA1, 2 (degrees)\n")
        if dim != 2:
            envfil.write(f"{p['Nbeta']}             ! Nbeta\n")
            envfil.write(f"{' '.join(map(str, p['beta']))} /  ! beta1, beta2 (degrees) bearing angle fan\n")
        if dim == 2:
            envfil.write(f"{p['deltas']} 5500.0 101.0 ! deltas, box depth, box range\n")
        else:
            envfil.write(f"{p['deltas']} 101.0 101.0 5500.0 ! deltas, box X, box Y, box depth\n")
        if it in {'R', 'C'}:
            envfil.write('\'MS\' 1.0 100.0 0, ! \'Width Curvature\' epsMultiplier rLoop ISINGL (ignored)\n')
            envfil.write('1 4 \'P\' ! Nimage iBeamWindow Component\n')
    ssp = p['ssp']
    if st in {'Q', 'H'}:
        with open('test/in/' + env_name + '.ssp', 'w') as sspfil:
            if st == 'Q':
                NPts, Nr, Rmin, Rmax = ssp['NPts'], ssp['Nr'], ssp['Rmin'], ssp['Rmax']
                def gen_ssp(r, z):
                    return 1500.0 + 50.0 * r / Nr + 50.0 * z / NPts
                sspfil.write(str(Nr) + '\n')
                sspfil.write(' '.join('{:.2f}'.format(Rmin+(Rmax-Rmin)*r/(Nr-1)) for r in range(Nr)) + '\n')
                for z in range(NPts):
                    sspfil.write(' '.join('{:.2f}'.format(gen_ssp(r, z)) for r in range(Nr)) + '\n')
            else:
                Nx, Ny, Nz = ssp['Nx'], ssp['Ny'], ssp['Nz']
                xymin, xymax, zmax = ssp['xymin'], ssp['xymax'], ssp['zmax']
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
    
def get_default_p(rt):
    p = {}
    p['NSx'] = 2
    p['Sx'] = [-20.0, 20.0]
    p['NSy'] = 2
    p['Sy'] = [-20.0, 20.0]
    p['NSz'] = 2
    p['Sz'] = [347.0, 682.0]
    if rt in {'E', 'A', 'a'}:
        p['NRz'] = 2
        p['Rz'] = [1135.8, 1145.8]
        p['NRr'] = 2
        p['Rr'] = [37.2, 37.21]
        p['Ntheta'] = 2
        p['theta'] = [43.2, 43.9]
    else:
        p['NRz'] = 51
        p['Rz'] = [0.0, 5000.0]
        p['NRr'] = 51
        p['Rr'] = [0.0, 100.0]
        p['Ntheta'] = 10
        p['theta'] = [0.0, 360.0]
    p['Nalpha'] = 21
    p['alpha'] = [-51.2, 51.2]
    p['Nbeta'] = 7
    p['beta'] = [0.0, 360.0]
    p['deltas'] = 1000.0
    p['ssp'] = {
        'NPts': 3,
        'Nr': 5,
        'Rmin': -102.0,
        'Rmax': 102.0,
        'Nx': 4,
        'Ny': 6,
        'Nz': 5,
        'xymin': -150.0,
        'xymax': 150.0,
        'zmax': 5000.0,
    }
    return p
    
def gen_all_it_st_combos(passtxt, failtxt, dim, rt):
    p = get_default_p(rt)
    if dim == 4 and rt not in {'E', 'A', 'a'}:
        p['theta'] = [1.0, 361.0]
    it_list = ['G'] if rt == 'R' else infl_types.keys()
    for it in it_list:
        for st in ssp_types.keys():
            subnames = [dims[dim], run_types[rt], infl_types[it], ssp_types[st]]
            env_name = 'gen_' + '_'.join(subnames)
            title = ', '.join(subnames)
            write_env_etc(dim, rt, it, st, p, env_name, title)
            (passtxt if should_work(dim, rt, it, st) else failtxt).write(env_name + '\n')

def gen_coverage_tests():
    for dim in [2, 3, 4]:
        def gen_coverage_tests_type(type, rt_list):
            tdname = type + dims[dim]
            print('\n\n{}:\n'.format(tdname))
            with open('gen_{}_pass.txt'.format(tdname), 'w') as passtxt, \
                open('gen_{}_fail.txt'.format(tdname), 'w') as failtxt:
                for rt in rt_list:
                    gen_all_it_st_combos(passtxt, failtxt, dim, rt)
        gen_coverage_tests_type('ray', ['R'])
        gen_coverage_tests_type('tl', ['C', 'S', 'I'])
        gen_coverage_tests_type('eigen', ['E'])
        gen_coverage_tests_type('arr', ['A', 'a'])

def gen_perf_ray_tests():
    for dim in [2, 3, 4]:
        env_name_base = 'genperf_ray' + dims[dim] + '_numrays'
        with open(env_name_base + '.txt', 'w') as txt:
            for Nalpha in [0x10, 0x40, 0x100, 0x400, 0x1000, 0x4000, 0x10000]:
                if dim != 2 and Nalpha > 1024:
                    break
                title = 'Increasing num rays ' + str(Nalpha)
                env_name = env_name_base + '_' + str(Nalpha)
                txt.write(env_name + '\n')
                rt, it, st = 'R', 'G', 'C'
                p = get_default_p(rt)
                p['Nalpha'] = Nalpha
                write_env_etc(dim, rt, it, st, p, env_name, title)
    #
    for dim in [2, 3, 4]:
        env_name_base = 'genperf_ray' + dims[dim] + '_long'
        Nalpha = 1000
        Nbeta = 1
        with open(env_name_base + '.txt', 'w') as txt:
            title = 'Long rays bc short max step'
            env_name = env_name_base + '_smallstep'
            txt.write(env_name + '\n')
            rt, it, st = 'R', 'G', 'C'
            p = get_default_p(rt)
            p['Nalpha'], p['Nbeta'] = Nalpha, Nbeta
            p['deltas'] = (5.0 if dim == 2 else 20.0)
            write_env_etc(dim, rt, it, st, p, env_name, title)
            #
            title = 'Long rays bc complex SSP'
            env_name = env_name_base + '_bigssp'
            txt.write(env_name + '\n')
            rt, it, st = 'R', 'G', ('Q' if dim == 2 else 'H')
            p = get_default_p(rt)
            p['Nalpha'], p['Nbeta'] = Nalpha, Nbeta
            p['ssp']['NPts'] = 200
            p['ssp']['Nr'] = 200
            p['ssp']['Nx'] = 200
            p['ssp']['Ny'] = 200
            p['ssp']['Nz'] = 200
            write_env_etc(dim, rt, it, st, p, env_name, title)

def gen_perf_tl_tests():
    for dim in [2, 3, 4]:
        env_name_base = 'genperf_tl' + dims[dim] + '_numrays'
        with open(env_name_base + '.txt', 'w') as txt:
            for Nalpha in [0x10, 0x40, 0x100, 0x400, 0x1000, 0x4000, 0x10000, 0x40000, 0x100000, 0x400000]:
                title = 'Increasing num rays ' + str(Nalpha)
                env_name = env_name_base + '_' + str(Nalpha)
                txt.write(env_name + '\n')
                rt, it, st = 'C', 'G', 'C'
                p = get_default_p(rt)
                p['NSx'] = p['NSy'] = p['NSz'] = 1
                p['Sx'] = p['Sy'] = [0.0]
                p['Sz'] = [823.4]
                p['NRz'] = p['NRr'] = 101
                p['Ntheta'] = 101
                p['theta'] = [0.0, 5.0]
                p['Nalpha'] = Nalpha
                p['alpha'] = [-20.1, 20.1]
                p['Nbeta'] = 1
                p['beta'] = [2.3]
                write_env_etc(dim, rt, it, st, p, env_name, title)
    #
    for dim in [2, 3, 4]:
        env_name_base = 'genperf_tl' + dims[dim] + '_numsources'
        with open(env_name_base + '.txt', 'w') as txt:
            for NS in range(1, 11):
                if dim != 2 and NS > 4:
                    break
                title = 'Increasing num sources '
                if dim == 2:
                    nsources = NS * NS
                    title += str(nsources)
                else:
                    title += str(NS) + 'x' + str(NS) + 'x' + str(NS)
                env_name = env_name_base + '_' + str(NS)
                txt.write(env_name + '\n')
                rt, it, st = 'C', 'G', 'C'
                p = get_default_p(rt)
                if dim == 2:
                    p['NSx'] = p['NSy'] = 1
                    p['NSz'] = nsources
                else:
                    p['NSx'] = p['NSy'] = p['NSz'] = NS
                p['Sz'] = [102.3, 4871.8]
                p['NRz'] = p['NRr'] = 101
                p['Ntheta'] = 101
                p['theta'] = [0.0, 5.0]
                p['Nalpha'] = 5000
                p['Nbeta'] = 1
                p['beta'] = [2.3]
                write_env_etc(dim, rt, it, st, p, env_name, title)
    #
    for Nalpha in [300, 20000]:
        for dim in [2, 3, 4]:
            env_name_base = 'genperf_tl' + dims[dim] + '_numreceivers'
            env_name_base += 'fewrays' if Nalpha == 300 else 'manyrays'
            with open(env_name_base + '.txt', 'w') as txt:
                for NR in [10, 30, 100, 300, 1000, 3000, 10000]:
                    if dim != 2 and NR > 300:
                        break
                    if dim == 2 and Nalpha == 20000 and NR == 10000:
                        break
                    title = 'Increasing num receivers with few rays, ' + str(NR) + 'x' + str(NR)
                    if dim != 2: title += 'x' + str(NR)
                    env_name = env_name_base + '_' + str(NR)
                    txt.write(env_name + '\n')
                    rt, it, st = 'C', 'G', 'C'
                    p = get_default_p(rt)
                    p['NSx'] = p['NSy'] = p['NSz'] = 1
                    p['Sx'] = p['Sy'] = [0.0]
                    p['Sz'] = [823.4]
                    p['NRz'] = p['NRr'] = p['Ntheta'] = NR
                    p['theta'] = [0.0, 30.0]
                    p['Nalpha'] = Nalpha
                    p['Nbeta'] = 1
                    p['beta'] = [14.2]
                    write_env_etc(dim, rt, it, st, p, env_name, title)

def gen_perf_arr_tests():
    for hits in [False, True]:
        for bigfield in [False, True]:
            for dim in [2, 3, 4]:
                env_name_base = 'genperf_arr' + dims[dim]
                env_name_base += '_hits' if hits else '_nohits'
                env_name_base += '_big' if bigfield else '_small'
                with open(env_name_base + '.txt', 'w') as txt:
                    for Nalpha in [0x10, 0x40, 0x100, 0x400, 0x1000, 0x4000, 0x10000, 0x40000, 0x100000, 0x400000]:
                        title = 'Arrivals, '
                        title += 'hits, ' if hits else 'no hits, '
                        title += 'big field, ' if bigfield else 'small field, '
                        title += str(Nalpha) + ' rays'
                        env_name = env_name_base + '_' + str(Nalpha)
                        txt.write(env_name + '\n')
                        rt, it, st = 'A', 'G', 'C'
                        p = get_default_p(rt)
                        p['NSx'] = p['NSy'] = p['NSz'] = 1
                        p['Sx'] = p['Sy'] = [0.0]
                        p['Sz'] = [823.4]
                        if bigfield:
                            p['Rz'] = [2000.0, 3000.0]
                            p['Rr'] = [31.2, 32.2]
                            p['theta'] = [38.1, 42.1]
                            # p['NRz'] = p['NRr'] = 3000 if dim == 2 else 100
                            # p['Ntheta'] = 1 if dim == 2 else 100
                            p['NRz'] = p['NRr'] = 300 if dim == 2 else 50
                            p['Ntheta'] = 1 if dim == 2 else 50
                        else:
                            p['Rz'] = [1135.8, 1145.8]
                            p['Rr'] = [37.2, 37.21]
                            p['theta'] = [43.2, 43.9]
                            p['NRz'] = p['NRr'] = 100 if dim == 2 else 20
                            p['Ntheta'] = 1 if dim == 2 else 20
                        Nbeta = 8 if dim != 2 else 1
                        p['Nalpha'] = Nalpha // Nbeta
                        p['Nbeta'] = Nbeta
                        p['alpha'] = [150.0, 170.0] if dim == 2 and not hits else [-51.2, 51.2]
                        p['beta'] = [30.0, 50.0] if hits else [182.3, 185.8]
                        write_env_etc(dim, rt, it, st, p, env_name, title)
                

gen_coverage_tests()
gen_perf_ray_tests()
gen_perf_tl_tests()
gen_perf_arr_tests()
