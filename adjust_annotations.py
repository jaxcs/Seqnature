#! /usr/bin/env python

# adjust_annotations.py:
# Given a set of strain-specific adjustment (offset) files, such as produced by
# the companion script build_reference_from_sanger_vcfs.py, adjust a gtf file
# to accommodate Indels.

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



import sys
import os
from optparse import OptionParser
import bisect

def parseOptions():
    usage = 'USAGE: %prog [options] annotation-file strain-identifier'

    parser = OptionParser(usage=usage, version="0.1")

    parser.add_option('-c', '--chr-column', dest='chrCol',
                      help='Column with chromosome info. (Required, one-based.)')
    parser.add_option('-C', '--comments-file', dest='comments',
                      help='File to write adjustment comments to. (Default: no comments file.)')
    parser.add_option('-d', '--directory', dest='dir', default='.',
                      help='path to directory containing offset files (default: current directory)')
    parser.add_option('-e', '--end-column', dest='endCol',
                      help='The column containing the position of the end of the feature. (Required, one-based.)')
    parser.add_option('-n', '--near', dest='near',
                      help='Report if a feature position is within near bases of an indel. (Optional, default is no report.)')
    parser.add_option('-o', '--output-file', dest='ofn',
                      help='Output filename (Optional, default: stdout)')
    parser.add_option('-p', '--position-columns', dest='posCols',
                      help='Comma-separated list of other columns with positions to be adjusted. (Optional, one-based.)')
    parser.add_option('-q', '--quiet', dest='quiet', action='store_true',
                      help='Operate quietly. Do not report chromosomes and contigs not in the Indel vcf file.')
    parser.add_option('-s', '--start-column', dest='startCol',
                      help='The column containing the position of the end of the feature. (Required, one-based.)')
    parser.add_option('-t', '--type-column', dest='typeCol',
                      help='Column with feature type info. (Optional, one-based.)')

    (opts, args) = parser.parse_args()

    errors = ''
    if len(args) != 2:
        errors += 'You must specify the input annotation file and strain!'

    # Sanity check our options and arguments.
    # Mandatory options are reference, indels and snps
    if not opts.chrCol:
        errors += '\n    You must specify a chromosome column.'
    if not opts.startCol:
        errors += '\n    You must specify the start position column.'
    if not opts.endCol:
        errors += '\n    You must specify the end position column.'
    strains = ['129P2', '129S1', '129S5', 'A_J', 'AKR', 'BALB', 'C3H', 'C57BL', 'CAST', 'CBA', 'DBA', 'LP_J', 'NOD', 'NZO', 'PWK', 'SPRET', 'WSB']
    if not args[1] in strains:
        errors += '\n    The strain you specified: ' + strain + ' is not a supported strain.'
    if errors:
        parser.error(errors)
    return (opts, args)

def readOffsets(dir, strain):
    chrs = []
    for n in range(1,20):
        chrs.append(str(n))
    chrs.append('X')
    
    chrOffs = {}
    chrKeys = {}
    for chr in chrs:
        fn = '%s_offsets_chr%s.txt' % (strain, chr)
        fp = os.path.join(dir, fn)
        try:
            f = open(fp)
            f.close()
        except IOError:
            print >> sys.stderr, 'Could not open', fp, '\nExiting...'
            sys.exit(1)
        offsets = {}
        keys = []
        for line in open (fp):
            line = line.rstrip()
            parts = line.split('\t')
            if len(parts) != 2:
                print >> sys.stderr, 'ERROR: file', fn, 'has bad line', line
                continue
            offsets[int(parts[0])] = int(parts[1])
            keys.append(int(parts[0]))
        chrOffs[chr] = offsets
        chrKeys[chr] = sorted(keys)

    return chrOffs, chrKeys


def adjustPosition(pos, offsets, keys):
    pos = int(pos)
    offsetIdx = bisect.bisect(keys, pos)

    if offsetIdx ==  0:
        # The pos is before the first offset entry, so the offset is 0.
        eventPos = 0
        offset = 0
    else:
        eventPos = keys[offsetIdx - 1]
        offset = offsets[eventPos]
    if offsetIdx < 2:
        offsetDelta = offset
    else:
        # Our offsetDelta is the difference between the current
        # offset and the previous one. If negative, this is a deletion.
        offsetDelta = offset - offsets[keys[offsetIdx - 2]]

    # Do some "quality" (exceptional events) checks.
    # Is the position in a deletion?
    inDeletion =  offsetDelta < 0 and pos < eventPos - offsetDelta
    

    # The new pos is pos + offset; the original (reference) pos is 'pos'.
    #if pos == 5193483 or pos == 5194069:
    #    print pos+ offset, pos, eventPos, inDeletion, offsetIdx, offset, offsetDelta
    #    if pos == 5194069:
    #        sys.exit()
    return pos + offset, pos, eventPos, inDeletion

def checkNear(chr, startRef, endRef, refPos, adjPos, eventPos, featureType, near, outcf):
    if near and abs(refPos - eventPos) <= near:
        if type(adjPos) == type(''):
            print >> sys.stderr, 'ERROR Bad adjPos:', adjPos
            sys.exit(1)
        #print >> sys.stderr, type(startRef), type(endRef), type(refPos), type(adjPos), type(near)
        if outcf:
            if featureType:
                print >> outcf, featureType,
            print >> outcf, '%s\t%d\t%d\t%d\t%d\tPosition near event (%d bases)' % (chr, startRef, endRef, refPos, adjPos, near)

alreadyCommented = set()
def doComment(feature, chr, start, end, ref, adj, format, comment, outcf):
    global alreadyCommented
    if not outcf:
        return
    if feature:
        key = feature
    key = key + chr + str(start) + str(end) + comment
    if key in alreadyCommented:
        return
    alreadyCommented.add(key)
    if feature:
        print >> outcf, feature,
    # Chr \t RefStart \t RefEnd \t RefLoc \t AdjLoc \t Description
    print >> outcf, format % (chr, start, end, ref, adj, comment)


otherChrs = []
lineNo = 0
def processLine(line, chrCol, startCol, endCol, posCols, chrOffs, chrKeys, typeCol, near, outf, outcf):
    global otherChrs, quiet
    line = line.rstrip()
    parts = line.split('\t')
    
    chr = parts[chrCol]
    if chr[:3] == 'chr':
        chr == chr[3:]
    try:
        offsets = chrOffs[chr]
    except KeyError:
        if chr != "###" and chr not in otherChrs and not quiet:
            otherChrs.append(chr)
            # Print diagnostic that we're ignoring chrs, once per chr.
            print >> sys.stdout, '# chromosome', chr, 'is not in the Sanger data; passing data through to the output file without adjusting the positions.'
        print >> outf, line
        return

    offset = None
    changeInFeature = False
    featureConsumed = False
    inDeletion = []
    nearHit = []
    keys = chrKeys[chr]

    # Process the start and end positions, then any extra positions.
    # for each perform some checks.
    startNewPos, startRefPos, startEventPos, startInDeletion = adjustPosition(parts[startCol], offsets, keys)
    

    endNewPos, endRefPos, endEventPos, endInDeletion = adjustPosition(parts[endCol], offsets, keys)

    if typeCol:
        featureType = parts[typeCol] + '\t'
    else:
        featureType = None

    if startInDeletion and endInDeletion and startEventPos == endEventPos:
        # The entire feature was deleted.
        # Set the start and end positions to zero.
        startNewPos = endNewPos = 0
        
        # and log it.
        doComment(featureType, chr, startRefPos, endRefPos, '', '', '%s\t%d\t%d\t%s\t%s\t%s', 'Feature Deleted', outcf)

    else:
        # At least part of the feature is in the strain. Check for both the
        # start and end.  And note that these are NOT mutually exclusive.
        # The start could be in one del and the end could be in another.
        # (Heaven help us if that really happens in biology...)
        if startInDeletion:
            # Adjust it to be the first non-deleted base. Add this segment's
            # offset to the event pos.
            startNewPos = startEventPos + (startNewPos - startRefPos)
            # and log it.
            comment = 'Start in deletion, adjusted to first non-deleted base'
            doComment(featureType, chr, startRefPos, endRefPos, startRefPos, startNewPos,
                      '%s\t%d\t%d\t%d\t%d\t%s', comment, outcf)
        if endInDeletion:
            # Adjust it to be the last base before the deletion.
            # Have to find out the offset between the strain and the reference just before this
            # event.
            testRefPos = endEventPos - 1;
            testNewPos, testRefPos, testEventPos, testInDeletion = adjustPosition(testRefPos, offsets, keys)
            # A quick sanity test.  Did we get a different offset?
            if (testNewPos - testRefPos) == (endNewPos - endRefPos):
                print >> sys.stderr, 'ERROR: We came up with the same offset!', chr, startRefPos, endRefPos
            endNewPos = testNewPos
            # and log it.
            comment = 'End in deletion, adjusted to last base before deletion'
            doComment(featureType, chr, startRefPos, endRefPos, endRefPos, endNewPos,
                      '%s\t%d\t%d\t%d\t%d\t%s', comment, outcf)

    if startEventPos != endEventPos:
        comment = 'Feature contains indel'
        doComment(featureType, chr, startRefPos, endRefPos, startRefPos, startNewPos,
                  '%s\t%d\t%d\t%d\t%d\t%s', comment, outcf)

    if endNewPos < startNewPos:
        print >> sys.stderr, 'ERROR: end before start', chr + ':' + str(endNewPos), 'in feature at', str(startNewPos) + '-' + str(endNewPos), '(reference coords:', str(startRefPos) + '-' + str(endRefPos) + ')'

    # Check start and end for "nearness"
    if near:
        checkNear(chr, startRefPos, endRefPos, startRefPos, startNewPos, startEventPos, featureType, near, outcf)
        checkNear(chr, startRefPos, endRefPos, endRefPos, endNewPos, endEventPos, featureType, near, outcf)
    
    for col in posCols:
        newPos, refPos, eventPos, inDeletion = adjustPosition(parts[col], offsets, keys)
        if inDeletion:
            comment = 'Position in deletion, adjusted to base after deletion'
            doComment(featureType, chr, startRefPos, endRefPos, refPos, newPos,
                      '%s\t%d\t%d\t%d\t%d\t%s', comment, outcf)
            newPos = eventPos + (newPos - refPos)
            parts[col] = str(newPos)

        if newPos < startNewPos:
            print >> sys.stderr, 'ERROR position before start', chr + ':' + str(newPos), 'in feature at', str(startNewPos) + '-' + str(endNewPos), '(reference coords:', str(startRefPos) + '-' + str(endRefPos) + ')'

        if newPos > endNewPos:
            print >> sys.stderr, 'ERROR position after end', chr + ':' + str(newPos), 'in feature at', str(startNewPos) + '-' + str(endNewPos), '(reference coords:', str(startRefPos) + '-' + str(endRefPos) + ')'

        # Is this position within (the user specified) 'near' bases of an
        # indel?
        if near:
            checkNear(chr, startRefPos, endRefPos, refPos, newPos, eventPos, featureType, near, outcf)

    # Update the file
    parts[startCol] = str(startNewPos)
    parts[endCol] = str(endNewPos)
    print >> outf, '\t'.join(parts)


def updateHeader(strain, outf):
    print >> outf, '# The coordinates in this file have been adjusted for strain %s' % strain
    print >> outf, "# based on Sanger's short Indel vcf file for that strain."
    print >> outf, '# This file is based on:\n#'


def emitCommentsHeader(typeCol, outcf):
    header = 'Chr\tRefStart\tRefEnd\tRefLoc\tAdjLoc\tDescription'
    if typeCol:
        header = 'FeatureType\t' + header
    print >> outcf, header

quiet = None

def main():
    global quiet
    (opts, args) = parseOptions()

    anns = args[0]
    strain = args[1]
    
    if opts.near:
        near = int(opts.near)
    else:
        near = None

    chrOffs, chrKeys   = readOffsets(opts.dir, strain)
    quiet = opts.quiet

    # Capture the interesting columns and make them 0-based. (Input
    # as 1-based.)
    chrCol = int(opts.chrCol) - 1
    typeCol = int(opts.typeCol) - 1
    startCol = int(opts.startCol) - 1
    endCol = int(opts.endCol) - 1

    if chrCol == -1:
        print >> sys.stderr, 'Columns are counted from 1, not zero.'
        sys.exit(1)
    if opts.posCols:
        posCols = opts.posCols.split(',')
    else:
        posCols = []
    for n in range(len(posCols)):
        posCols[n] = int(posCols[n]) - 1

    try:
        f = open(anns)
        f.close()
    except IOError:
        print >> sys.stderr, 'Could not open', anns, '\nExiting...'
        sys.exit(1)

    outf = sys.stdout  # Main annotations output
    outcf = None
    if opts.ofn:
        try:
            outf = open(opts.ofn, 'w')
        except IOError:
            print >> sys.stderr, 'Could not open', opts.ofn, 'for output.\nExiting...'
    if opts.comments:
        try:
            #Try to open the comments file.
            outcf = open(opts.comments, 'w')
        except IOError:
            print >> sys.stderr, 'Could not open comments file', opts.comments, 'for output. Continuing processing without comments.'
            outcf = None

    headerProcessed = False
    headerUpdated = False
    if outcf:
        emitCommentsHeader(opts.typeCol, outcf)

    for line in open(anns):
        if not headerProcessed:
            if line[0] == '#':
                print >> outf, line,
                if not headerUpdated:
                    updateHeader(strain, outf)
                    headerUpdated = True
                continue
            else:
                headerProcessed = True

        processLine(line, chrCol, startCol, endCol, posCols, chrOffs, chrKeys, typeCol, near, outf, outcf)

    if outf != sys.stdout:
        outf.close()
    if outcf:
        outcf.close()

# Do it!
main()
