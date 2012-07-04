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

class ExonSet(object):
    def __init__(self, exons):
        self.exons = exons  # a list of exons

    def search(self, cur_ex):
        """
        assumes that cur_ex will always be found
        in self.exons
        """
        l = 0
        r = len(self.exons) - 1
        while l < r:
            m = int((l + r) / 2)
            if self.exons[m] == cur_ex:
                return self.exons[m]
            if self.exons[m] > cur_ex:
                r = m - 1
            else:
                l = m + 1
        return self.exons[r]        

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

class GeneLocus(ExonSet):
    def __init__(self, exs):
        super(GeneLocus, self).__init__(exs)
        self.transcript_ids = None  # dictionary of each transcript's exons

def build_gene_loci(tr_exs):
    """
    tr_exs has to be in the format of the output of 
        gtf_0.get_transcripts_exons
    """
    chrs = {}
    
    exon_t_ids = {}  
    # each entry is transcript_id: set of 
    # transcript's names that contain this exon 
    
    for tr_name, exs in tr_exs.iteritems():
        chrm = chrs.setdefault(exs[0].seqname, set([]))
        for ex in exs:
            cur_ex = Exon(ex.start, ex.end)
            chrm.add(cur_ex)
            exon_t_ids.setdefault(cur_ex, set([])).add(tr_name)
            
    for ex, tr_names in exon_t_ids.iteritems():
        exon_t_ids[ex] = list(tr_names)
            
    for chr_name, chrm in chrs.iteritems():
        cur_chrom = ExonSet(sorted(list(chrm), cmp=Exon.exon_cmp))
        chrs[chr_name] = cur_chrom
    
    mod_tr_exs = {}
    
    class TmpTr(ExonSet):
        def __init__(self, exs, chr_name):
            super(TmpTr, self).__init__(exs)
            self.chr_name = chr_name
    
    for tr_name, exs in tr_exs.iteritems():
        mod_exs = []
        for ex in exs:
            mod_exs.append(Exon(ex.start, ex.end, connect=False))
        mod_exs = sorted(mod_exs, cmp=Exon.exon_cmp)
        mod_tr_exs[tr_name] = TmpTr(mod_exs, exs[0].seqname)
    
    for tr_name, tmp_tr in mod_tr_exs.iteritems():
        chr_exs = []
        for tmp_ex in tmp_tr.exons:
            chr_exs.append(chrs[tmp_tr.chr_name].search(tmp_ex))
        if len(chr_exs) > 1:
            chr_exs[0].right_exons.append(chr_exs[1])
            chr_exs[-1].left_exons.append(chr_exs[-2])
        for chr_ex, i in izip(chr_exs[1:-1], xrange(1, len(chr_exs) - 1)):
            chr_ex.left_exons.append(chr_exs[i - 1])
            chr_ex.right_exons.append(chr_exs[i + 1])
    
    gene_loci = {}
    for chr_name, chrm in chrs.iteritems():
        cur_gene_loci = []
        locusized_exs = set([])
        for ex in chrm.exons:
            if ex not in locusized_exs:
                
                def get_gene_locus(cur_ex):
                    my_locus = [cur_ex]
                    locusized_exs.add(cur_ex)
                    for l_ex in cur_ex.left_exons:
                        if l_ex not in locusized_exs:
                            my_locus.extend(get_gene_locus(l_ex))
                    for r_ex in cur_ex.right_exons:
                        if r_ex not in locusized_exs:
                            my_locus.extend(get_gene_locus(r_ex))
                    return my_locus
                
                cur_gene_locus = get_gene_locus(ex)
                if cur_gene_locus:
                    cur_gene_locus = sorted(cur_gene_locus, cmp=Exon.exon_cmp)
                    cur_gl = GeneLocus(cur_gene_locus)
                    cur_ts = set([])
                    for ex in cur_gene_locus:
                        for t_id in exon_t_ids[ex]:
                            cur_ts.add(t_id)
                    t_exs = {}
                    for t_id in cur_ts:
                        cur_t_exs = []
                        for ex in mod_tr_exs[t_id].exons:
                            cur_t_exs.append(cur_gl.search(ex))
                        t_exs[t_id] = sorted(cur_t_exs, cmp=Exon.exon_cmp)
                    cur_gl.transcript_ids = t_exs 
                    cur_gene_loci.append(cur_gl)
        gene_loci[chr_name] = cur_gene_loci
    return gene_loci


