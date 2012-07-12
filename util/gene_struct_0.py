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

import sys

class ExonSet(object):
    def __init__(self, exons):
        self.exons = exons  # a list of exons

class Exon(object):
    def __init__(self, start, end, connect=True):
        # python indexing
        self.start = start
        self.end = end
        if connect:
            self.left_exons = set([])
            self.right_exons = set([])
    
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
    
    def __repr__(self):
        return "start: %d, end: %d" % (self.start, self.end)
    
    def __str__(self):
        return self.__repr__()
 
class GeneLocus(ExonSet):
    def __init__(self, exs):
        super(GeneLocus, self).__init__(exs)
        self.t_exs = {}  # dictionary of each transcript's exons
    
    def count_chains(self):
        dfsed = set([])
        
        def count_chain_DFS(cur_ex):
            if cur_ex in dfsed:
                return
            dfsed.add(cur_ex)
        
        count_chain_DFS(self.exons[0])
                

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
    
    max_ex_num = 0
    
    tmp_chrm_exs = {}
    for chrm_name, exs in chrm_exs.iteritems():
        exs = sorted(list(exs), cmp=Exon.exon_cmp)
        fixed_exs = []
        i = 0
        cur_block = set([])
        while i < len(exs):
            cur_ex = Exon(exs[i].start, exs[i].end)
            cur_block.add(exs[i].start)
            cur_block.add(exs[i].end)
            i += 1
            while i < len(exs) and cur_ex.overlap(exs[i]):
                cur_block.add(exs[i].start)
                cur_block.add(exs[i].end)
                if exs[i].end > cur_ex.end:
                    cur_ex.end = exs[i].end
                i += 1
            cur_block = sorted(list(cur_block))
            for j in xrange(0, len(cur_block) - 1):
                fixed_exs.append(Exon(cur_block[j], cur_block[j + 1]))
            cur_block = set([])
        tmp_chrm_exs[chrm_name] = fixed_exs
        if len(fixed_exs) > max_ex_num:
            max_ex_num = len(fixed_exs)
    
    chrm_exs = tmp_chrm_exs
    
    sys.setrecursionlimit(2 * max_ex_num)
    
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
    
    ex_t_ids = {}
    
    for t_name, fixed_exs in t_fixed_exs.iteritems():
        for ex in fixed_exs:
            ex_t_ids.setdefault((ex, tr_exs[t_name][0].seqname),
                                set([])).add(t_name)
        if len(fixed_exs) > 1:
            fixed_exs[0].right_exons.add(fixed_exs[1])
            fixed_exs[1].left_exons.add(fixed_exs[0])
            fixed_exs[-1].left_exons.add(fixed_exs[-2])
            fixed_exs[-2].right_exons.add(fixed_exs[-1])
        for i in xrange(1, len(fixed_exs) - 1):
            fixed_exs[i].right_exons.add(fixed_exs[i + 1])
            fixed_exs[i + 1].left_exons.add(fixed_exs[i])
            fixed_exs[i].left_exons.add(fixed_exs[i - 1])
            fixed_exs[i - 1].right_exons.add(fixed_exs[i])
    
    gene_loci = {}
    
    for chrm_name, exs in chrm_exs.iteritems():
        found_exs = set([])
        cur_gloci = []
        for ex in exs:
            if ex not in found_exs:
                
                def find_gene_locus_exons(cur_ex):
                    found_exs.add(cur_ex)
                    cur_exs = [cur_ex]
                    for l_ex in cur_ex.left_exons:
                        if l_ex not in found_exs:
                            cur_exs.extend(find_gene_locus_exons(l_ex))
                    for r_ex in cur_ex.right_exons:
                        if r_ex not in found_exs:
                            cur_exs.extend(find_gene_locus_exons(r_ex))
                    return cur_exs
                
                gl_exs = sorted(find_gene_locus_exons(ex), Exon.exon_cmp)
                
                cur_ts = set([])
                for ex1 in gl_exs:
                    cur_ts = cur_ts | ex_t_ids[(ex1, chrm_name)]
                cur_gl = GeneLocus(gl_exs)
                for t_name in cur_ts:
                    cur_gl.t_exs[t_name] = t_fixed_exs[t_name]
                    
                cur_gloci.append(cur_gl)
        gene_loci[chrm_name] = cur_gloci
        
    """
    return a dictionary
        chrm_name: list of sorted gene_loci
    """
    return gene_loci



