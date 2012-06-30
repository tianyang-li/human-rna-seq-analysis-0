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

from utils.blat_0 import blat_reader

def main():
    blat_file1, blat_file2 = None, None
    try:
        opts, args = getopt.getopt(sys.argv[1:], '',
                                   ['blat1=', 'blat2='])
    except getopt.GetoptError as err:
        print >> sys.stderr, str(err)
        sys.exit(1)
    for opt, arg in opts:
        if opt == '--blat1':
            blat_file1 = arg
        if opt == '--blat2':
            blat_file2 = arg
    if (not blat_file1 
        or not blat_file2):
        print >> sys.stderr, "missing"
        sys.exit(1)
    
        
    
if __name__ == '__main__':
    main()
