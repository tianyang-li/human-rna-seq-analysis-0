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

import pysam

from util.gtf_0 import get_transcripts_exons
from util.gene_struct_0 import build_gene_loci

def main():
    gtf_file = None
    bam_file = None
    try:
        opts, _ = getopt.getopt(sys.argv[1:], '', ['gtf=', 'bam='])
    except getopt.GetoptError as err:
        print >> sys.stderr, str(err)
        sys.exit(1)
    for opt, arg in opts:
        if opt == '--gtf':
            gtf_file = arg
        if opt == '--bam':
            bam_file = arg
    if (not gtf_file 
        or not bam_file):
        print >> sys.stderr, "misisng"
        sys.exit(1)
    
    tr_exs = get_transcripts_exons(gtf_file)
    gene_loci = build_gene_loci(tr_exs)
    
    #TODO: analyze contigs structure

if __name__ == '__main__':
    main()

