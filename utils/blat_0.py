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

class BlatEntry(object):
    """
    0 matches - Number of bases that match that aren't repeats
    1 misMatches - Number of bases that don't match
    2 repMatches - Number of bases that match but are part of repeats
    3 nCount - Number of 'N' bases
    4 qNumInsert - Number of inserts in query
    5 qBaseInsert - Number of bases inserted in query
    6 tNumInsert - Number of inserts in target
    7 tBaseInsert - Number of bases inserted in target
    8 strand - '+' or '-' for query strand. For translated alignments, second '+'or '-' is for genomic strand
    9 qName - Query sequence name
    10 qSize - Query sequence size
    11 qStart - Alignment start position in query
    12 qEnd - Alignment end position in query
    13 tName - Target sequence name
    14 tSize - Target sequence size
    15 tStart - Alignment start position in target
    16 tEnd - Alignment end position in target
    17 blockCount - Number of blocks in the alignment (a block contains no gaps)
    18 blockSizes - Comma-separated list of sizes of each block
    19 qStarts - Comma-separated list of starting positions of each block in query
    20 tStarts - Comma-separated list of starting positions of each block in target
    """
    def __init__(self, line):
        line = line.strip()
        entries = line.split("\t")
        self.matches = int(entries[0])
        self.mis_matches = int(entries[1])
        self.rep_matches = int(entries[2])
        self.N_count = int(entries[3])
        self.Q_num_insert = int(entries[4])
        self.Q_base_insert = int(entries[5])
        self.T_num_insert = int(entries[6])
        self.T_base_insert = int(entries[7])
        self.strand = entries[8]
        self.Q_name = entries[9]
        self.Q_size = int(entries[10])
        self.Q_start = int(entries[11])
        self.Q_end = int(entries[12])
        self.T_name = entries[13]
        self.T_size = int(entries[14])
        self.T_start = int(entries[15])
        self.T_end = int(entries[16])
        self.block_count = int(entries[17])
        self.block_sizes = entries[18].split(",")
        self.block_sizes = map(lambda s: int(s), self.block_sizes)
        self.Q_starts = entries[19].split(",")
        self.Q_starts = map(lambda s:int(s), self.Q_starts)
        self.T_starts = entries[20].split(",")
        self.T_starts = map(lambda s:int(s), self.T_starts)

def blat_reader(blat_file):
    with open(blat_file, 'r') as fin:
        for line in fin:
            yield BlatEntry(line)


