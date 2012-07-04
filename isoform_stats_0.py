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
    
    for chr_gl in gene_loci.itervalues():
        for gl in chr_gl:
            gl_ex_num = len(gl.exons)
            left_end_diff = 0
            right_end_diff = 0
            both_ends_diff = 0
            for tr_exs in gl.transcript_ids.itervalues():
                left_diff = False
                right_diff = False
                if tr_exs[0].start != gl.exons[0].start:
                    left_diff = True
                    left_end_diff += 1 
                if tr_exs[-1].end != gl.exons[-1].end:
                    right_diff = True
                    right_end_diff += 1
                if left_diff and right_diff:
                    both_ends_diff += 1
            print gl_ex_num, left_end_diff, right_end_diff, both_ends_diff

if __name__ == '__main__':
    main()

