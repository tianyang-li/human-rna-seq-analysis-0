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

class Chr(object):
    def __init__(self):
        self.exons = []

class Exon(object):
    def __init__(self, start, end):
        self.start = start
        self.end = end
        self.left_exons = []
        self.right_exons = []
    
    def get_len(self):
        self.end - self.start
    
    def __cmp__(self, other):
        if self.start == other.start:
            return self.end - other.end
        return self.start - other.start
    
    def __hash__(self):
        return hash((self.start, self.end))
    
    @staticmethod
    def cmp(a, b):
        return a.__cmp__(b)
    

def build_gene_loci(tr_exs):
    """
    tr_exs has to be in the format of the output of 
        gtf_0.get_transcripts_exons
    """
    chrs = {}
    for exs in tr_exs.itervalues():
        chr = chrs.setdefault(exs[0].seqname, set([]))
        for ex in exs:
            chr.add(Exon(ex.start, ex.end))
    for chr in chrs.itervalues():
        chr = sorted(list(chr), cmp=Exon.cmp)
        
    
    
