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

"""
Structure is as GFF, so the fields are: 
    < seqname > 
    < source > 
    < feature > 
    < start > 
    < end > 
    < score > 
    < strand > 
    < frame > 
    [attributes] 
    [comments]
"""

class GTFEntry(object):
    def __init__(self, line):
        pos_hash = line.find("#")
        if pos_hash == -1:
            self.comments = ""
            entries = line.split("\t")
        else:
            self.comments = line[pos_hash + 1:]
            entries = line[:pos_hash].split("\t")
        self.seqname = entries[0]
        self.source = entries[1]
        self.feature = entries[2]
        self.start = int(entries[3])
        self.end = int(entries[4])
        if entries[5] == ".":
            self.score = None
        else:
            self.score = float(entries[5])
        self.strand = entries[6]
        if entries[7] == ".":
            self.frame = None
        else:
            self.frame = int(entries[7])
        self.attributes = entries[8]
    
    def __str__(self):
        if not self.score:
            score = "."
        else:
            score = "%f" % self.score
        if self.frame == None:
            frame = "."
        else:
            frame = "%d" % self.frame
        line = ("%s\t%s\t%s\t"
               "%d\t%d\t%s\t"
               "%s\t%s\t%s") % (self.seqname,
                      self.source,
                      self.feature,
                      self.start,
                      self.end,
                      score,
                      self.strand,
                      frame,
                      self.attributes)
        if not self.comments:
            return line
        line = "%s#%s" % (line, self.comments)
        return line
 
def gtf_reader(gtf_file):
    with open(gtf_file, 'r') as fin:
        for line in fin:
            yield GTFEntry(line.strip())



