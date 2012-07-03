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

from itertools import izip
import sys

from util.exception_0 import StrException

class Chr(object):
    def __init__(self, exons):
        self.exons = exons
    
    def search(self, cur_ex):
        l = 0
        r = len(self.exons)
        while l < r:
            m = int((l + r) / 2)
            if self.exons[m] == cur_ex:
                return self.exons[m]
            if self.exons[m] > cur_ex:
                r = m - 1
            else:
                l = m + 1
        if l > r:
            raise StrException("l > r in Chr.search")
        if self.exons[m] != cur_ex:
            raise StrException("self.exons[m] != cur_ex in Chr.search")
        return self.exons[m]

class Exon(object):
    def __init__(self, start, end, connect=True):
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
    for chr_name, chr in chrs.iteritems():
        chrs[chr_name] = Chr(sorted(list(chr), cmp=Exon.exon_cmp))
    
    mod_tr_exs = {}
    
    class TmpTr(object):
        def __init__(self, exs, chr_name):
            self.exs = exs
            self.chr_name = chr_name
    
    for tr_name, exs in tr_exs.iteritems():
        mod_exs = []
        for ex in exs:
            mod_exs.append(Exon(ex.start, ex.end, connect=False))
        mod_exs = sorted(mod_exs, cmp=Exon.exon_cmp)
        mod_tr_exs[tr_name] = TmpTr(mod_exs, exs[0].seqname)
    
    for tr_name, tmp_tr in mod_tr_exs.iteritems():
        chr_exs = []
        for tmp_ex in tmp_tr.exs:
            try:
                chr_exs.append(chrs[tmp_tr.chr_name].search(tmp_ex))
            except StrException as err:
                print >> sys.stderr, str(err)
                return
        if len(chr_exs) > 1:
            chr_exs[0].right_exons.append(chr_exs[1])
            chr_exs[-1].left_exons.append(chr_exs[-2])
        for chr_ex, i in izip(chr_exs[1:-1], xrange(1, len(chr_exs) - 1)):
            chr_ex.left_contigs = chr_exs[i - 1]
            chr_ex.right_contigs = chr_exs[i + 1]
    
