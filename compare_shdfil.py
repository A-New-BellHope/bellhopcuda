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
from math import isnan, isfinite

thresh_ulp_print = 4000
thresh_ulp_alarm = 100000
thresh_abs_print = 1e-7
thresh_abs_alarm = 1e-5
thresh_rel_print = 1e-4
thresh_rel_alarm = 1e-2

def compare_files(cxxf, forf):
    cxxdata = cxxf.read()
    fordata = forf.read()
    if len(cxxdata) != len(fordata):
        print('File lengths differ: {} vs {}'.format(len(cxxdata), len(fordata)))
        sys.exit(1)
    
    reclen = 0
    def read_rec_int(data, rec, val):
        start = 4*(rec*reclen+val)
        return struct.unpack('I', data[start:start+4])[0]
    cxxreclen = read_rec_int(cxxdata, 0, 0)
    forreclen = read_rec_int(fordata, 0, 0)
    if cxxreclen != forreclen:
        print('Record lengths differ: {} vs {}'.format(cxxreclen, forreclen))
        sys.exit(1)
    reclen = cxxreclen
    if len(cxxdata) % reclen != 0:
        print('File length {} is not a multiple of the record length {}'.format(
            len(cxxdata), reclen))
        sys.exit(1)
    
    Nfreq = read_rec_int(fordata, 2, 0)
    Ntheta = read_rec_int(fordata, 2, 1)
    NSx = read_rec_int(fordata, 2, 2)
    NSy = read_rec_int(fordata, 2, 3)
    NSz = read_rec_int(fordata, 2, 4)
    NRz = read_rec_int(fordata, 2, 5)
    NRr = read_rec_int(fordata, 2, 6)
    PlotType = read_rec_int(fordata, 1, 0)
    assert Nfreq == read_rec_int(cxxdata, 2, 0)
    assert Ntheta == read_rec_int(cxxdata, 2, 1)
    assert NSx == read_rec_int(cxxdata, 2, 2)
    assert NSy == read_rec_int(cxxdata, 2, 3)
    assert NSz == read_rec_int(cxxdata, 2, 4)
    assert NRz == read_rec_int(cxxdata, 2, 5)
    assert NRr == read_rec_int(cxxdata, 2, 6)
    assert PlotType == read_rec_int(cxxdata, 1, 0)
    assert NRr * 2 <= reclen
    isTL = (PlotType & 0xFFFF) == 0x4C54 # 'TL' (only write first and last Sx/Sy)
    irre = (PlotType == 0x65727269) #'irre' (gular)
    assert isTL or irre or (PlotType == 0x74636572) #'rect' (ilin)
    rcvrgridsz = Ntheta * (1 if irre else NRz) * reclen # reclen is normally NRr (*2 for complex)
    filesz = NSx * NSy * NSz * rcvrgridsz * 4 + 4 * 10 * reclen
    if len(cxxdata) != filesz:
        print('NSx {} NSy {} NSz {} Ntheta {} NRz {} NRr {} Nfreq {}'.format(
            NSx, NSy, NSz, Ntheta, NRz, NRr, Nfreq))
        print('reclen {} rcvrgridsz {} irregular {}\nPred filesz {} actual {}'.format(
            reclen, rcvrgridsz, irre, filesz, len(cxxdata)))
        print('Invalid file size')
        sys.exit(1)

    print('FOR title: {}'.format(struct.unpack('80s', fordata[4:84])[0].decode('ascii')))
    print('CXX title: {}'.format(struct.unpack('80s', cxxdata[4:84])[0].decode('ascii')))
    
    for rec in range(3, 10):
        N = {3: Nfreq, 4: Ntheta, 5: 2 if isTL else NSx, 6: 2 if isTL else NSy,
            7: NSz, 8: NRz, 9: NRr}[rec]
        type, l = ('d', 8) if rec == 3 else ('f', 4)
        for i in range(N):
            start = l*(rec*reclen+i)
            forf = struct.unpack(type, fordata[start:start+l])[0]
            cxxf = struct.unpack(type, cxxdata[start:start+l])[0]
            if isnan(forf) or isnan(cxxf) or (forf != cxxf and (cxxf - forf) / forf > 1e-8):
                print('Pos data did not match: rec {} element {}: FOR {} CXX {}'.format(
                    rec, i, forf, cxxf))
                sys.exit(1)
        for i in range((N*l//2), reclen):
            fori = read_rec_int(fordata, rec, i)
            cxxi = read_rec_int(cxxdata, rec, i)
            if not(fori == cxxi == 0):
                print('Zeroes after pos data did not match: rec {} element {} FOR {:08X} CXX {:08X}'.format(
                    rec, i, fori, cxxi))
                sys.exit(1)

    errors = 0
    maxerrors = 100
    for rec in range(10, len(cxxdata) // reclen):
        for Rr in range(0, NRr):
            addr = rec * reclen * 4 + Rr * 8
            cxxd = cxxdata[addr:addr+8]
            ford = fordata[addr:addr+8]
            if cxxd == ford:
                continue
            #
            d = rec - 10
            if irre:
                Rz = 0
            else:
                Rz = d % NRz
                d //= NRz
            Sz = d % NSz
            d //= NSz
            theta = d % Ntheta
            d //= Ntheta
            Sy = d % NSy
            d //= NSy
            Sx = d % NSx
            d //= NSx
            if d != 0:
                raise RuntimeError('Indexing failed')
            p_addr = 'Src XYZ {:3},{:3},{:3} Rcvr ThZR {:3},{:3},{:3} {:08X}: '.format(
                Sx, Sy, Sz, theta, Rz, Rr, addr)
            #
            cir, cii = struct.unpack('II', cxxd)
            fir, fii = struct.unpack('II', ford)
            cfr, cfi = struct.unpack('ff', cxxd)
            ffr, ffi = struct.unpack('ff', ford)
            fmt_flt = '{:12.6}'
            fmt_flt2 = '({},{}) / ({},{}) | '.format(fmt_flt, fmt_flt, fmt_flt, fmt_flt)
            p_flt = fmt_flt2.format(cfr, cfi, ffr, ffi)
            if isnan(cfr) != isnan(ffr) or isnan(cfi) != isnan(ffi):
                print('\n\n=============== NAN DETECTED ===============\n{}{}'.format(p_addr, p_flt))
                sys.exit(1)
            ulpr, ulpi = abs(cir - fir), abs(cii - fii)
            absr, absi = abs(cfr - ffr), abs(cfi - ffi)
            relr = absr / abs(ffr) if ffr != 0 else absr
            reli = absi / abs(ffi) if ffi != 0 else absi
            #
            p_int = '{:08X} {:08X} / {:08X} {:08X} | '.format(cir, cii, fir, fii)
            p_ulp = 'ULP ({:4},{:4}) | '.format(ulpr, ulpi)
            p_abs = 'ABS ({},{}) | '.format(fmt_flt, fmt_flt).format(absr, absi)
            p_rel = 'REL ({:9.6}%,{:9.6}%) | '.format(relr*100.0, reli*100.0)
            # If the absolute error is miniscule, it doesn't matter what the relative
            # or ULP error are
            if absr < 5e-9:
                ulpr = 0
                relr = 0.0
            if absi < 5e-9:
                ulpi = 0
                reli = 0.0
            if ulpr > thresh_ulp_alarm or ulpi > thresh_ulp_alarm or \
                absr > thresh_abs_alarm or absi > thresh_abs_alarm or \
                relr > thresh_rel_alarm or reli > thresh_rel_alarm:
                if errors == 0:
                    print('\nERROR: Extremely large error(s) detected:')
                if errors < maxerrors:
                    print(p_addr + p_flt + p_abs + p_rel)
                errors += 1
            if errors > 0:
                continue
            #
            errs = ''
            show_int, show_flt = False, False
            if ulpr > thresh_ulp_print or ulpi > thresh_ulp_print:
                errs += p_ulp
                show_int = True
            if absr > thresh_abs_print or absi > thresh_abs_print:
                errs += p_abs
                show_flt = True
            if relr > thresh_rel_print or reli > thresh_rel_print:
                errs += p_rel
                show_flt = True
            if not show_int and not show_flt:
                continue
            if show_flt:
                errs = p_flt + errs
            if show_int:
                errs = p_int + errs
            print(p_addr + errs)
        
        for addr in range(rec * reclen + NRr * 8, (rec+1) * reclen, 4):
            cxxd = cxxdata[addr:addr+4]
            ford = fordata[addr:addr+4]
            if cxxd != b'\0\0\0\0' or ford != b'\0\0\0\0':
                print('Files contain non-zero padding at {:08X}'.format(addr))
                sys.exit(1)

    if errors >= maxerrors:
        print('\nand {} more extremely large error(s)'.format(errors - maxerrors))
    if errors > 0:
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) not in {2, 3}:
        print('Usage: python3 compare_shdfil.py MunkB_Coh [cxx1/cxxmulti/cuda]')
        print('No paths, no .shd')
        sys.exit(1)

    shdfil = sys.argv[1]
    if len(sys.argv) == 3:
        comparisons = [sys.argv[2]]
    else:
        comparisons = ['cxx1', 'cxxmulti', 'cuda']

    for c in comparisons:
        with open('test/FORTRAN/{}.shd'.format(shdfil), 'rb') as forf:
            cxxfile = 'test/{}/{}.shd'.format(c, shdfil)
            try:
                with open(cxxfile, 'rb') as cxxf:
                    print('TL / shade file comparison FORTRAN vs. {}:'.format(c))
                    compare_files(cxxf, forf)
            except FileNotFoundError:
                print('{} not found, skipping {}'.format(cxxfile, c))
