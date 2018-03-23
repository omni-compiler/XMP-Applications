#!/usr/bin/env python

import math
import sys

STEP = 100

FIELD = 5  # total energ
#FIELD = 8  # pressure
#FIELD = 6  # temperature

EPS = 1.0e-8


def usage():
    print 'usage: %s master_file file' % sys.argv[0]


def read_val(file):
    f = open(file, 'r')
    for line in f:
        field = line.split()
        if field[0][0] == '#': continue
        if int(field[0]) == STEP: return float(field[FIELD])


if __name__ == '__main__':

    if len(sys.argv) != 3:
        usage()
        sys.exit(1)

    print
    print 'Checking result ...'

    val0 = read_val(sys.argv[1])
    val  = read_val(sys.argv[2])
    diff = math.fabs((val - val0) / val0)

    if diff < EPS:
        print 'OK (relative error = %e)' % diff
        sys.exit(0)
    else:
        print 'NG (relative error = %e)' % diff
        sys.exit(1)

    
