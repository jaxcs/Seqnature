#! /usr/bin/env python

# reconstructNCBI.py
# Re-order a fa file, e.g., NCBIM37.fa, into chromosome order, and split
# off Y, MT and the unplaced NT_* contigs into separate files.  This is 
# a companion program to build_individualized_genome.py.  Having the 
# reference ordered in this way greatly speeds up creation of each
# animal's custom reference.

#   Author: Al Simons 
#   Copyright (C) 2012  The Jackson Laboratory
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU Affero General Public License as
#   published by the Free Software Foundation, either version 3 of the
#   License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU Affero General Public License for more details.
#
#   You should have received a copy of the GNU Affero General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.

usage="""
 USAGE: reconstruct_NCBI.py reference-file.fa

 Output is written to three files in the current working directory.
 The file names are based on the input .fa name.  Assume the
 input file is named NCBIM37.fa.  The output files would be
   - NCBIM37_ordered.fa
   - NCBIM37_contigs.fa
   - NCBIM37_Y.fa

 The program for which this is a preprocessor, build_individualized_genome.py,
 takes in the basename, optionally with a path.  In our example, this 
 would be NCBIM37 or /some/path/to/NCBIM37

 About chromosome names: This program assumes NCBI-style chromosome
 names "1", "2", "X", etc., not mm9-style names, e.g., "chr1".
"""

import sys, os


if len(sys.argv) != 2:
    print >> sys.stderr, usage
    sys.exit(1)

# Define the order in which we want the chrs...
chrs = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
for n in range(len(chrs)):
    chrs[n] = str(chrs[n])
chrs.append('X')
prefix = os.path.splitext(os.path.basename(sys.argv[1]))[0]

def contigFilename():
    global prefix
    return prefix+'_contigs.fa'

def chrFilename(chr):
    global prefix
    return prefix + '_' + chr + '.fa'

of = None
contigFile = open(contigFilename(), 'w')
for line in open(sys.argv[1]):
    if line[0] == '>':
        if of and of != contigFile:
            of.close()
        if line[2] == 'T':
            # line[2] == 'T' gets both >NT_* and >MT
            of = contigFile
        else:
            chr = line.split()[0][1:]
            of = open(chrFilename(chr), 'w')
    of.write(line)

of.close()
contigFile.close()

# Note that we have created a Y chromosome file, and we are NOT going to
# include that into the merged, ordered file.  Also we don't merge in
# the contigs.

of = open(prefix+'_ordered.fa', 'w')
for chr in chrs:
    for line in open(chrFilename(chr)):
        print >> of, line,
    os.remove(chrFilename(chr))
