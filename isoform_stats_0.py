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

from util.gtf_0 import get_transcripts_exons
from util.gene_struct_0 import build_gene_loci

def main():
    gtf_file = None
    try:
        opts, _ = getopt.getopt(sys.argv[1:], '', ['gtf='])
    except getopt.GetoptError as err:
        print >> sys.stderr, str(err)
        sys.exit(1)
    for opt, arg in opts:
        if opt == '--gtf':
            gtf_file = arg
    if not gtf_file:
        print >> sys.stderr, "misisng"
        sys.exit(1)
    
    tr_exs = get_transcripts_exons(gtf_file)
    gene_loci = build_gene_loci(tr_exs)
    
    for chr_name, chr_gl in gene_loci.iteritems():
        for gl in chr_gl:
            chain_tot = len(gl.ex_chains)
            for exs in gl.t_exs.itervalues():
                chain_set = set([])
                for ex in exs:
                    chain_set.add(gl.ex2chain[ex])
                left_no_distinct = False
                right_no_distinct = False
                if exs[0].start != gl.ex2chain[exs[0]].exons[0].start:
                    left_no_distinct = True
                if exs[-1].end != gl.ex2chain[exs[-1]].exons[-1].end:
                    right_no_distinct = True
                print chain_tot, len(chain_set),
                if left_no_distinct:
                    print 1,
                else:
                    print 0,
                if right_no_distinct:
                    print 1,
                else:
                    print 0,
                if left_no_distinct and right_no_distinct:
                    print 1,
                else:
                    print 0,
                print

if __name__ == '__main__':
    main()

