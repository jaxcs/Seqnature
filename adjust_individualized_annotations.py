#! /usr/bin/env python

# adjust_DO_annotations.py:

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


# Given a set of DO animal-specific adjustment (offset) files, such as 
# produced by two runs (A and B) of the companion script
# build_DO_from_sanger_vcfs.py, create an animal specific gtf file as follows:
# 1. Chromosomes will designated <name>{L,R}.
# 2. Genes will have the code letters of the strains present in the overall
#    range of the gene, appended to the gene name after an underscore.
# 3. Transcripts will have the code letters of the strains contributing to
#    the transcripts exons, appended to the transcript name.

# To do this we also read in the animal's strain extents (probabilities) file.

# Since neither genes nor transcripts are explicitly recorded with ranges in
# the gtf file, we walk the gtf file twice.  Once to determine the strains
# contributing to each gene/transcript, and once to write the file with
# updated names and positions.  We do it in this order because the strain
# extent are recorded in reference coordinates.


import sys
import os
from optparse import OptionParser
import bisect

# En/Dis-able some logging to stdout / stderr:
DEBUG=False

def debug(msg):
    global DEBUG
    if DEBUG:
        print >> sys.stderr, msg

def parseOptions():
    usage = 'USAGE: %prog [options] annotation-file DO-identifier L-or-R'

    parser = OptionParser(usage=usage, version="0.1")

    parser.add_option('-c', '--chr-column', dest='chrCol',
                      help='Column with chromosome info. (Required, one-based.)')
    parser.add_option('-C', '--comments-file', dest='comments',
                      help='File to write adjustment comments to. (Default: no comments file.)')    
    parser.add_option('-d', '--dir', dest='dir', default='.',
                      help='path to directory containing offset files (default: current directory)')
    parser.add_option('-e', '--end-column', dest='endCol', default='5', 
                      help='The column containing the position of the end of the feature. (Optional, one-based. Default: 5, the GTF standard column.)')
    parser.add_option('-n', '--near', dest='near',
                      help='Report if a feature position is within near bases of an indel. (Optional, default is no report.)')
    parser.add_option('-o', '--output-file', dest='ofn',
                      help='Output filename (Optional, default: stdout)')
    parser.add_option('-p', '--position-columns', dest='posCols',
                      help='Comma-separated list of other columns with positions to be adjusted. (Optional, one-based.)')
    parser.add_option('-s', '--start-column', dest='startCol', default='4',
                      help='The column containing the position of the end of the feature. (Optional, one-based. Default: 4, the GTF standard column.)')
    parser.add_option('-t', '--type-column', dest='typeCol',
                      help='Column with feature type info. (Optional, one-based.)')
    parser.add_option('-x', '--extents', dest='extents',
                      help='File containing the strain extents.')
    (opts, args) = parser.parse_args()

    errors = ''
    if len(args) != 3:
        errors += 'You must specify the input annotation file, DO id and A/B designator!'

    # Sanity check our options and arguments.
    # Mandatory options are reference, indels and snps
    if not opts.chrCol:
        errors += '\n    You must specify a chromosome column.'
    if not opts.startCol:
        errors += '\n    You must specify the start position column.'
    if not opts.endCol:
        errors += '\n    You must specify the end position column.'
    if not opts.extents:
        errors += '\n    You must specify the file with strain extents.'

    if errors:
        parser.error(errors)
    return (opts, args)


def readStrainBounds(fN):
    chrs = []
    for n in range(1,20):
        chrs.append(str(n))
    chrs.append('X')
    
    chrStrains = {}
    chrKeys = {}

    headerSkipped = False
    for line in open(fN):
        if not headerSkipped:
            headerSkipped = True
            continue
        (chr, offset, strain) = line.split(',')
        chr = chr.strip(' \t')
        offset = int(offset.strip(' \t'))
        strain = strain.strip(' \t\r\n')
        if chr not in chrStrains:
            chrStrains[chr] = {}
            chrKeys[chr] = []
        chrStrains[chr][offset] = strain
        chrKeys[chr].append(offset)

    for chr in chrs:
        chrKeys[chr] = sorted(chrKeys[chr])

    return chrStrains, chrKeys


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
            debug('Could not open' + str(fp) + 'assuming all B6.')
            continue
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

def getRangeStrains(start, end, strains, keys):
    # Return the codes for all the strains that contribute to this range.
    # For exons, this will almost always be a single strain; for genes,
    # it may well be two or even more.
    startIdx = bisect.bisect(keys, start)
    endIdx = bisect.bisect(keys, end)

    if startIdx == 0 or endIdx == 0:
        # There's been a problem.  Should be no entry before the first.
        print >> sys.stderr, "ERROR looking up strain info at positions", start, end
        sys.exit(1)

    retStrains = []
    for n in range(startIdx-1, endIdx):
        transitionPos = keys[n]
        retStrains.append(strains[transitionPos])
    return retStrains

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
    return pos + offset, pos, eventPos, inDeletion

def checkNear(chr, startRef, endRef, refPos, adjPos, eventPos, featureType, near, outcf):
    if near and abs(refPos - eventPos) <= near:
        if type(adjPos) == type(''):
            print >> sys.stderr, 'ERROR Bad adjPos:', adjPos
            sys.exit(1)
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


female = False

def processLine(line, chrCol, startCol, endCol, posCols, chrOffs, chrKeys, typeCol, near, outf, outcf, geneNames, transcriptNames, proteinNames, lr):
    global female
    line = line.rstrip()
    parts = line.split('\t')
    
    chr = parts[chrCol]
    if chr[:3] == 'chr':
        chr == chr[3:]
    if female and  'Y' in chr:
        return
    try:
        offsets = chrOffs[chr]
    except KeyError:
        # If a chr is entirely B6 it will not have an offsets file,
        # and won't appear in the chrOffs list.
        # The correct thing to do is simply pass all the lines for that
        # chr through unmodified except for properly munging the chr name,
        # and gene and transcript ids.
        
        # The wrinkle (another??) is that NT_, MT and Y also don't have 
        # an offsets file and should be passed through completely unchanged.
        if chr[0] == 'Y' or chr[0] == 'N' or chr[0] == 'M':
            print >> outf, line
            return
        parts[chrCol] = chr + lr
        notes = parts[8].split(';')
        for n in range(len(notes)):
            for id in ["gene_id", "transcript_id"]:
                if id in notes[n]:
                    # The strain identifier for B6 is "B", hard coded in this line...
                    notes[n] = ' %s "%sB%s"' % (id, notes[n].split(' ')[2][1:-1], lr)
        parts[8] = ';'.join(notes)        
        print >> outf, '\t'.join(parts)
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

    # Fix up the "notes" column (assumed to be last). The gene_id is the
    # first part, and the transcript_id is the second.  Both are updated.
    notes = parts[-1].split(';')
    if len(notes) < 2:
        print >> sys.stderr, 'ERROR: notes are too short!', parts
        sys.exit(1)
    if notes[0][:8] != ' gene_id':
        print >> sys.stderr, 'ERROR: gene_id not first in notes...', notes
        sys.exit(1)
    if notes[1][:14] != ' transcript_id':
        print >> sys.stderr, 'ERROR: transcript_id not second in notes...', notes
        sys.exit(1)

    g = notes[0].split('"')[1]
    t = notes[1].split('"')[1]
    g = geneNames[g]
    t = transcriptNames[t]
    notes[0] = ' gene_id "%s"' % (g)
    notes[1] = 'transcript_id "%s"' % (t)

    p = notes[-1]
    if p[:11] == ' protein_id':
        p = p.split('"')[1]
        notes[-1] = 'protein_id "%s"' % (proteinNames[p])

    parts[-1] = '; '.join(notes)
        
    
        
    # Update the file
    parts[chrCol] = chr + lr
    parts[startCol] = str(startNewPos)
    parts[endCol] = str(endNewPos)
    print >> outf, '\t'.join(parts)

class NameAnnotator:
    def __init__(self, name, chr, lr):
        self.chr = chr
        self.name = name
        self.start = 999999999999  # Only used for genes
        self.end = 0   # Only used for genes
        self.strains = [lr]
    def addStrains(self, strains):
        for strain in strains:
            if strain not in self.strains:
                self.strains.append(strain)
    def expandRange(self, start, end):
        if start < self.start:
            self.start = start
        if end > self.end:
            self.end = end
    def __str__(self):
        return self.name + ''.join(sorted(self.strains))

def annotateGenesTranscriptsProteins(anns, chrCol, typeCol, startCol, endCol, chrStrains, chrStrainKeys, lr):
    genes = {}
    transcripts = {}
    proteins = {}

    notesCol = -1
    for line in open(anns):
        line = line.rstrip()
        parts = line.split('\t')
        # Only exons contribute to the transcripts and define the extent
        # of the genes...
        if parts[typeCol] != 'exon':
            continue
        chr = parts[chrCol]
        start = int(parts[startCol])
        end = int(parts[endCol])
        notes = parts[notesCol].split(';')
        gene_id = notes[0].split('"')[1]
        transcript_id = notes[1].split('"')[1]
        protein_id = notes[-1]
        if protein_id[:11] == ' protein_id':
            protein_id = protein_id.split('"')[1]
        else:
            protein_id = None

        if chr in chrStrains:
            strains = getRangeStrains(start, end, chrStrains[chr], chrStrainKeys[chr])
        else:
            strains = ''
        # Annotate the gene names,
        if gene_id not in genes:
            genes[gene_id] = NameAnnotator(gene_id, chr, lr)
        gene = genes[gene_id]
        gene.expandRange(start, end)
        # Since we want to annorate a gene with all the strains it spans,
        # even those not in its exons, we'll add the strains later.

        # Annotate transcripts
        if transcript_id not in transcripts:
            transcripts[transcript_id] = NameAnnotator(transcript_id, chr, lr)
        transcript = transcripts[transcript_id]
        transcript.addStrains(strains)

        # And proteins
        if protein_id:
            if protein_id not in proteins:
                proteins[protein_id] = NameAnnotator(protein_id, chr, lr)
            protein = proteins[protein_id]
            protein.addStrains(strains)
        
    # We have now processed all the exons in the gtf.  We can now
    # 1) Update the strains in all the genes, and 
    # 2) Update the genes and transcripts hashes with the annotated names,
    #    throwing away the annotators in the process.
    geneKeys = genes.keys()
    for gene in geneKeys:
        g = genes[gene]
        c = g.chr
        if c in chrStrains:
            geneStrains = getRangeStrains(g.start, g.end, chrStrains[c], chrStrainKeys[c])
            g.addStrains(geneStrains)
        genes[gene] = str(g)

    # Clean up the transcripts
    transcriptKeys = transcripts.keys()
    for transcript in transcriptKeys:
        transcripts[transcript] = str(transcripts[transcript])

    # and the proteins
    proteinKeys = proteins.keys()
    for protein in proteinKeys:
        proteins[protein] = str(proteins[protein])

    return genes, transcripts, proteins

def updateHeader(strain, outf):
    print >> outf, '# The coordinates in this fle have been adjusted for DO animal genotype identified as', strain
    print >> outf, "# based on Sanger's short Indel vcf files for the DO founder strains."
        
    print >> outf, '# This file is based on:\n#'


def emitCommentsHeader(typeCol, outcf):
    if not outcf:
        return
    header = 'Chr\tRefStart\tRefEnd\tRefLoc\tAdjLoc\tDescription'
    if typeCol:
        header = 'FeatureType\t' + header
    print >> outcf, header


def main():
    global female
    (opts, args) = parseOptions()

    anns = args[0]
    lr = args[2]
    animalId = args[1] + '_' + lr
    
    if animalId[0] == 'F':
        female = True

    if opts.near:
        near = int(opts.near)
    else:
        near = None

    chrOffs, chrKeys   = readOffsets(opts.dir, animalId)
    chrStrains, chrStrainKeys = readStrainBounds(opts.extents)

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

    geneNames, transcriptNames, proteinNames = annotateGenesTranscriptsProteins(anns, chrCol, typeCol, startCol, endCol, chrStrains, chrStrainKeys, lr)

    outf = sys.stdout  # Main annotations output
    outcf = None # Comments about the annotation translation process.
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
            print >> sys.stderr, 'Could not open comments file', opts.comments, 'for output. Continuing without comments.'
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

        processLine(line, chrCol, startCol, endCol, posCols, chrOffs, chrKeys, typeCol, near, outf, outcf, geneNames, transcriptNames, proteinNames, lr)

    if outf != sys.stdout:
        outf.close()
    if outcf:
        outcf.close()

# Do it!
main()

