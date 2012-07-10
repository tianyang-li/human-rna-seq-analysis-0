#!/usr/bin/env python

#  Copyright (C) 2012 Tianyang Li
#  tmy1018@gmail.com
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License

import getopt
import sys
from itertools import izip

from Bio import SeqIO

def main():
    fin1, fin2 = None, None
    fmt = None
    out_prefix = None
    try:
        opts, _ = getopt.getopt(sys.argv[1:], 'f:1:2:', ['prefix='])
    except getopt.GetoptError as err:
        print >> sys.stderr, str(err)
        sys.exit(1)
    for opt, arg in opts:
        if opt == '-f':
            fmt = arg
        if opt == '-1':
            fin1 = arg
        if opt == '-2':
            fin2 = arg
        if opt == '--prefix':
            out_prefix = arg
            
    if (not fin1
        or not fin2
        or not fmt
        or not out_prefix):
        print >> sys.stderr, "missing"
        sys.exit(1)
        
    fout1 = open("%s_1.%s" % (out_prefix, fmt), 'w')
    fout2 = open("%s_2.%s" % (out_prefix, fmt), 'w')
    
    for rec1, rec2 in izip(SeqIO.parse(fin1, fmt), SeqIO.parse(fin2, fmt)):
        tmp_seq = "%s%s" % (str(rec1.seq), str(rec2.seq))
        if "N" not in tmp_seq and "n" not in tmp_seq:
            fout1.write(rec1.format(fmt))
            fout2.write(rec2.format(fmt))
    
    fout1.close()
    fout2.close()

if __name__ == '__main__':
    main()    


