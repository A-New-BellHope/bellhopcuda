import sys, struct

if len(sys.argv) != 2:
    print('Usage: python3 compare_shdfil.py MunkB_Coh')
    print('No path, no .shd')
    sys.exit(1)

thresh_ulp_print = 1000
thresh_ulp_alarm = 10000
thresh_abs_print = 1e-6
thresh_abs_alarm = 1e-4
thresh_rel_print = 5e-6
thresh_rel_alarm = 1e-4

cxxfile = 'temp/cxx/{}.shd'.format(sys.argv[1])
forfile = 'temp/FORTRAN/{}.shd'.format(sys.argv[1])

with open(cxxfile, 'rb') as cxxf, open(forfile, 'rb') as forf:
    cxxdata = cxxf.read()
    fordata = forf.read()
    if len(cxxdata) != len(fordata):
        print('File lengths differ: {} vs {}'.format(len(cxxdata), len(fordata)))
        sys.exit(0)
    
    cxxreclen = struct.unpack('I', cxxdata[0:4])[0]
    forreclen = struct.unpack('I', fordata[0:4])[0]
    if cxxreclen != forreclen:
        print('Record lengths differ: {} vs {}'.format(cxxreclen, forreclen))
        sys.exit(0)
    
    for w in range(10*cxxreclen):
        cxxd = cxxdata[w*4:(w+1)*4]
        ford = fordata[w*4:(w+1)*4]
        if cxxd != ford:
            print('{:08X}: CXX {:08X}  FOR {:08X}'.format(w*4,
                struct.unpack('I', cxxd)[0], struct.unpack('I', ford)[0]))
    
    for addr in range(4*10*cxxreclen, len(cxxdata), 8):
        cxxd = cxxdata[addr:addr+8]
        ford = fordata[addr:addr+8]
        if cxxd == ford:
            continue
        cir, cii = struct.unpack('II', cxxd)
        fir, fii = struct.unpack('II', ford)
        cfr, cfi = struct.unpack('ff', cxxd)
        ffr, ffi = struct.unpack('ff', ford)
        ulpr, ulpi = abs(cir - fir), abs(cii - fii)
        absr, absi = abs(cfr - ffr), abs(cfi - ffi)
        relr, reli = absr / ffr, absi / ffi
        p_addr = '{:08X}: '.format(addr)
        p_int = '{:08X} {:08X} / {:08X} {:08X} | '.format(cir, cii, fir, fii)
        p_ulp = 'ULP ({:4},{:4}) | '.format(ulpr, ulpi)
        fmt_flt = '{:12.6}'
        fmt_flt2 = '({},{}) / ({},{}) | '.format(fmt_flt, fmt_flt, fmt_flt, fmt_flt)
        p_flt = fmt_flt2.format(cfr, cfi, ffr, ffi)
        p_abs = 'ABS ({},{}) | '.format(fmt_flt, fmt_flt).format(absr, absi)
        p_rel = 'REL ({},{}) | '.format(fmt_flt, fmt_flt).format(relr, reli)
        if ulpr > thresh_ulp_alarm or ulpi > thresh_ulp_alarm or \
            absr > thresh_abs_alarm or absi > thresh_abs_alarm or \
            relr > thresh_rel_alarm or reli > thresh_rel_alarm:
            print('\nERROR: Extremely large error detected:\n' 
                + p_addr + p_int + p_flt + p_ulp + p_abs + p_rel + '\n')
            sys.exit(1)
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
