import os
import string

def getTotal(fname):
    rLine = ''
    if os.path.isfile(fname):
        f = open(fname, 'r')
        for line in f.readlines():
            key = line.split(' ')
            if key[0] == 'Total':
                d = key[3].split('d')
                rLine = d[0]+d[1]
        f.close()
    return rLine

of = open('.\\GTotal.csv', 'w')
for fn in range(0,99):
    fcode = str(fn).zfill(2)
    line = getTotal('.\\G' + fcode + '.log')
    if line != '':
        of.write(fcode + ',' + line)
of.close()


