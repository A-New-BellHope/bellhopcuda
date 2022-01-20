import sys, struct

if sys.version_info.major < 3:
    print('This is a python3 script')
    sys.exit(-1)

with open(sys.argv[1], 'rb') as f1, open(sys.argv[2], 'rb') as f2:
    d1 = f1.read()
    d2 = f2.read()

assert len(d1) == len(d2)

rleval = -1
rlelen = 0
np = 0
for a in range(len(d1) // 4):
    f1 = struct.unpack('<I', d1[4*a:4*a+4])[0]
    f2 = struct.unpack('<I', d2[4*a:4*a+4])[0]
    delta = f2 - f1
    if delta != rleval:
        if rleval >= 0:
            pass#print('{}x {}'.format(rlelen, rleval))
        rleval = delta
        rlelen = 1
    else:
        rlelen += 1
    if delta > 999 or delta < -99:
        if np != 0:
            print('')
        print('{}'.format(delta))
        np = 0
    else:
        print('{:3d} '.format(delta), end='')
        np += 1
    if np == 25:
        print('')
        np = 0

#print('{}x {}'.format(rlelen, rleval))
print('')
