import os
import string
import sys
inFName = sys.argv[1]
outFName = inFName[:-4]+'_G'
f = open(inFName, 'r')
line = f.readline()
keywords = line.split(',')

for line in f:
    v = line.split(',')
    of = open(outFName + v[0] + '.h', 'w')
    #of.write(keywords(2), values(2))
    i = 1
    lastK = ''
    s = ''
    for k in keywords[1:]:
        if k == 'End':
            break
        if not (k in ('index', 'type', 'value', 'dir', 'Phase')):
            if all(c in "0123456789.+-" for c in v[i]):
                value = v[i]
            else:
                value = '"' + v[i] + '"'
            if k == lastK:
                s += ' ' + value
            else:
                s += '\n' + k + ' ' + value
            lastK = k
        else:
            if k == 'index':
                s += '\n' + v[i]
            else:
                s += ' ' + '"' + v[i] + '"'
        i += 1
    of.write(s)
    of.close()
