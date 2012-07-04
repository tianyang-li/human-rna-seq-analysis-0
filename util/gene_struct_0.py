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

from __future__ import division

class ExonSet(object):
    def __init__(self, exons):
        self.exons = exons  # a list of exons

class Exon(object):
    def __init__(self, start, end, connect=True):
        # python indexing
        self.start = start
        self.end = end
        if connect:
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
    def exon_cmp(a, b):
        return a.__cmp__(b)
    
    def overlap(self, other):
        if self.end <= other.start:
            return False
        if self.start >= other.end:
            return False
        return True
    
    @staticmethod
    def exon_overlap(a, b):
        return a.overlap(b)

class GeneLocus(ExonSet):
    def __init__(self, exs):
        super(GeneLocus, self).__init__(exs)
        self.transcripts = {}  # dictionary of each transcript's exons

def build_gene_loci(tr_exs):
    """
    tr_exs has to be in the format of the output of 
        gtf_0.get_transcripts_exons
    """
    
    chrm_exs = {}
    
    for tr_ex in tr_exs.itervalues():
        chrm_name = tr_ex[0].seqname
        for ex in tr_ex:
            chrm_exs.setdefault(chrm_name, set([])).add(Exon(ex.start, ex.end))
    
    for chrm_name, exs in chrm_exs.iteritems():
        exs = sorted(list(exs), cmp=Exon.exon_cmp)
        fixed_exs = []
        i = 0
        cur_block = set([])
        for ex in exs:
            cur_block.add(ex.start)
            cur_block.add(ex.end)
        cur_block = sorted(list(cur_block))
        for i in xrange(0, len(cur_block) - 1):
            fixed_exs.append(Exon(cur_block[i], cur_block[i + 1]))
        chrm_exs[chrm_name] = fixed_exs
    
    t_fixed_exs = {}
    
    for tr_name, tr_ex in tr_exs.iteritems():
        
        def gtf_exon_cmp(a, b):
            if a.start == b.start:
                return a.end - b.end
            return a.start - b.start
        
        tr_ex = sorted(tr_ex, cmp=gtf_exon_cmp)
        chrm_name = tr_ex[0].seqname
        fixed_exs = chrm_exs[chrm_name]
        tr_start = tr_ex[0].start
        tr_end = tr_ex[-1].end
        
        def get_exon_by_start():
            l = 0
            r = len(fixed_exs) - 1
            while l < r:
                m = int((l + r) / 2)
                if fixed_exs[m].start == tr_start:
                    return m
                if fixed_exs[m].start > tr_start:
                    r = m - 1
                else:
                    l = m + 1
            return l
        
        start_exon = get_exon_by_start()
            
        def get_exon_by_end():
            l = 0
            r = len(fixed_exs) - 1
            while l < r:
                m = int((l + r) / 2)
                if fixed_exs[m].end == tr_end:
                    return m
                if fixed_exs[m].end > tr_end:
                    r = m - 1
                else:
                    l = m + 1
            return l
        
        end_exon = get_exon_by_end()
        
        t_fixed_exs[tr_name] = fixed_exs[start_exon:end_exon + 1]
    
    for fixed_exs in t_fixed_exs.itervalues():
        if len(fixed_exs) > 1:
            fixed_exs[0].right_exons.append(fixed_exs[1])
            fixed_exs[-1].left_exons.append(fixed_exs[-2])
        for i in xrange(1, len(fixed_exs) - 1):
            fixed_exs[i].right_exons.append(fixed_exs[i + 1])
            fixed_exs[i].left_exons.append(fixed_exs[i - 1])
    
    gene_loci = {}

    return gene_loci


