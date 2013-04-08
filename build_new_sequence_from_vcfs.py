#! /usr/bin/env python
#
# build_reference_from_sanger_vcfs.py
#
# Using SNP and indel VCF files from Sanger, produce reference genomes for
# seven of the eight CC founder strains.
#

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

def parseOptions():
    usage = 'USAGE: %prog [options] strain\n    Strain must be a column header in the two VCF files.'

    parser = OptionParser(usage=usage, version="1.0")

    parser.add_option('-d', '--output-dir', dest='dir', default='.',
                      help='Directory for output files (must already exist). (Default: current working directory.)')
    parser.add_option('-i', '--indels', dest='indels',
                      help='path to indels vcf file')
    parser.add_option('-o', '--output', dest='output',
                      help='File into which the genome is written. (Default: stdout)')
    parser.add_option('-r', '--reference', dest='reference',
                      help='Path to the reference file used as the base.')
    parser.add_option('-s', '--snps', dest='snps',
                      help='Path to the snps vcf file')

    (opts, args) = parser.parse_args()

    if len(args) != 1:
        print >> sys.stderr, args
        parser.error('You must specify the strain!')

    # Sanity check our options and arguments.
    # Mandatory options are reference, indels and snps
    errors = ''
    if not opts.indels:
        errors += '\n    You must specify an indels vcf file.'
    if not opts.reference:
        errors += '\n    You must specify a reference sequence file.'
    if not opts.snps:
        errors += '\n    You must specify a snps vcf file.'

    if errors:
        parser.error(errors)
    return (opts, args)

class Reference:
    # Even though it is not coded as one this is in fact a singleton.
    def __init__(self, chrs, inF):
        self.chrs = chrs
        # FIXME DEBUG
        print >> sys.stderr, chrs
        self.nextChr = ''
        self.currChr = ''
        self.chrLine = ''
        self.chrSeq = ''
        self.file=open(inF)
        self.inRewindSearch = False



    # FIXME DEBUG
    aksdebug = None
    def printChrs(self, msg):
        if not Reference.aksdebug:
            Reference.aksdebug = open('AKS_DEBUG.txt', 'w')
        print >> Reference.aksdebug, msg, self.chrs
        print >> sys.stderr, msg, self.chrs



    def findChr(self, chr):
        # Set up to read a particular chr
        while chr != self.nextChr:
            line = self.file.readline()
            if len(line) == 0:
                # Reached EOF.
                # If we're not already in a rewind search, we need to try;
                # we may have simply asked for things out of order.
                # If, however, we are in a rewind search, then the
                # requested chr isn't in the reference.  We need to bomb
                # the run. (sob... wipe away tear...)
                if not self.inRewindSearch:
                    self.inRewindSearch = True
                    self.file.seek(0)
                    self.findChr(chr)
                else:
                    print >> sys.stderr, 'Could not find chr', chr, 'in the reference. Exiting...'
                    sys.exit(1)

            if line[:1] == '>':
                # We've found a chromosome line.
                # These come in multiple forms. mm9 has:
                # >chr1
                # NCBI has:
                # >1 dna:chromosome chromosome...
                # We'll isolate the chromosome name, and stash the original
                # line so we can output it.
                #
                line = line.rstrip()
                thisChr = line[1:]
                if thisChr[:3] == 'chr':
                    thisChr = thisChr[3:]
                thisChr = thisChr.split(' ')[0]
                self.chrLine = '>' + str(thisChr)
                if thisChr == chr:
                    self.nextChr = thisChr
                    self.inRewindSearch = False

    def readChr(self, chr):
        seq = []
        if chr != self.nextChr:
            self.findChr(chr)
        self.currChr = self.nextChr
        while True:
            line = self.file.readline()
            if len(line) == 0:
                # EOF.  Save what we have...
                self.chrSeq = ''.join(seq)
                self.nextChr = 'EOF'
                return
            line = line.rstrip()
            if line[0] == '>':
                # Chromosome break
                if self.currChr == '':
                    # First time.
                    self.currChr = line[4:]
                else:
                    self.nextChr = line[4:]
                    self.chrSeq = ''.join(seq)
                    return
            else:
                # A normal sequence line
                seq.append(line)

    def getRange(self, chr, beg, end):
        # Retreives a string representing the reference from base pair
        # 'beg', through base pair 'end', inclusive.
        # Notice that pairs are one based, and the string is zero based.
        # that's why we use beg-1 and end, rather than the expected
        # beg and end+1.
        #
        if chr != self.currChr:
            if chr in self.chrs:
                # Indicate that we've seen this chr, by removing it
                # from the chrs list.  Then at the end, we'll output 
                # without modifications all those chrs we didn't see
                # while processing the SNPs and indels.
                del(self.chrs[self.chrs.index(chr)])

                # FIXME DEBUG
                self.printChrs("Removing chr" + chr)

            self.readChr(chr)
        return self.chrSeq[beg-1:end]

# End of class Reference...

# A common class for getting info from snps and indels.
class Result:
    def __init__(self, chr, pos, ref, allele, change, nextRefPos, isIndel):
        self.chr = chr
        self.pos = pos
        self.ref = ref
        self.allele = allele
        self.change = change
        self.nextRefPos = nextRefPos
        self.isIndel = isIndel

    def __str__(self):
        return ('I: ' if self.isIndel else 'S: ') + self.chr + ":" + str(self.pos) + ' ' + self.ref + '->' + self.allele + '(' + str(self.change) + ') nrp:' + str(self.nextRefPos)


#
# A class to track indels.  Use toResult() to get results.
# Members:
#     chr - self explanatory
#     position - This is NOT the position contained in the VCF file, which
#                was the last common base between the reference and the
#                allele.  This is some number of bases later, the position
#                of the first different base.  Consider the case of 
#                chr1 3000742, ref: GT, alleles G,GTT,GGTTTT.
#                The first allele is a one base deletion.  The second is
#                a one base insertion, 2 bases after the listed position.
#                the third is a change of the second base and then an
#                insertion of 4 bases. Position for this allele would
#                be one after the pos in the file.
#     allele -   The bases (for an insert) starting at the first different
#                base. Null string for a deletion.
#     change -   the effect of this allele on the new coordinates.
#                len(allele) - len(ref)
#     nextRefPos - The next base position that is to be taken from the
#                reference.
#
class Indel:
    def __init__(self, file, column):
        # Get the next line in the file with a call for this strain.
        # Return a list of the tab-sep line 
        first = True
        while first or parts[column] == '.':
            line = file.readline()
            if len(line) == 0:
                # We've reached EOF.
                self.chr = 'EOF'
                return 
            line = line.rstrip()
            parts = line.split('\t')
            first = False
        # Now we have a line with an indel for this strain.
        self.chr = parts[0]
        position = int(parts[1])
        ref = parts[3]
        alleleGroup = parts[4]
        call = parts[column]

        alleles = alleleGroup.split(',')
        parts1 = call.split(':')
        parts2 = parts1[0].split('/')
        if parts2[0] != parts2[1]:
            print >> sys.stderr, 'het call found at chr' + chr, 'position', position, call, 'Using first allele.'
        #Calls use 1-based allele references.
        try:
            self.allele = alleles[int(parts2[0])-1]
        except IndexError:
            print >> sys.stderr, 'Index error processing', parts
            sys.exit(1)
            

        #
        # The way the Sanger VCFs are created, the first base in both the
        # reference and the allele will be the same. Further bases _may_
        # be the same. Make position be the first different base and fix up the
        # allele to not contain the common bases.
        allele = self.allele
        minLen = min(len(ref), len(allele))
        diffPos = 0
        for n in range(minLen):
            if ref[n] != allele[n]:
                diffPos = n
                break

        # If diffPos == 0, the reference and allele matched throughout
        # their length.  This is either an insertion after the reference,
        # or a deletion.
        # Either way, the difference position is the base after the shorter
        # of the two.
        if diffPos == 0:
            diffPos = minLen
        self.position = position + diffPos
        self.allele = self.allele[diffPos:]
        self.ref = ref[diffPos:]

        # A somewhat common pattern is for a few bases at the beginning to be 
        # deleted, followed by many bases that are the same.  
        # This is necessitated by the multi-strain nature of the file. 
        # Some other strain may have deleted many more bases here.
        lenDiff = len(self.ref) - len(self.allele)
        if lenDiff > 0 and self.ref[lenDiff:] == self.allele:
            self.ref = self.ref[:lenDiff]
            self.allele = ''
        
        # Also, the change in position will be len(allele) - len(ref).  This
        # will be negative for a deletion, and positive for an insertion.
        self.change = len(self.allele) - len(self.ref)
        self.nextRefPos = self.position + len(self.ref)

    def toResult(self):
        return Result(self.chr, self.position, self.ref, self.allele, self.change, 
                      self.nextRefPos, True)

    def __str__(self):
        return self.chr + ':' + str(self.position) + ' ' + self.ref + '->' + self.allele + ' ' + str(self.change) + ' ' + str(self.nextRefPos)

    def __lt__(self, snp):
        if isinstance(snp, Snp):
            if self.chr != snp.chr:
                # Currently, the VCF files SEEM TO contain chrs 1-N and X, Y, in that
                # order.  THIS PROGRAM ASSUMES THAT'S TRUE...  Hope so.
                try:
                    return int(self.chr) < int(snp.chr)
                except ValueError:
                    # One or both is an alphabetic chr name.
                    # Figure out whether either is numeric. It is "less"
                    try:
                        int(self.chr)
                        # If we get here, then we are less. (self is numeric,
                        # so snp has to be alphabetic, which we regard as
                        # greater.)
                        return True
                    except:
                        pass
                    try:
                        int(snp.chr)
                        # And if we get here, self is greater.
                        return False
                    except:
                        Pass
                    # Both were alphabetic.  Do simple string sort.
                    return self.chr < snp.chr

            # Same chr, return based on position, already ints.
            return self.position < snp.position
        print >> sys.stderr, 'Must be called with a Snp argument.'
        return NotImplemented

#
# A class to track SNPs. 
# Members:
#     chr - self explanatory
#     position - the position of the snp from the file.
#     allele - the snp
#     change - 0 (snps are always a single base replacement; no length
#              differences.) (pseudo member; materialized when needed.)
#     nextRefPos - position + 1 (pseudo member)
class Snp:
    def __init__(self, file, column):
        # Get the next line in the file with a call for this strain.
        # Return a list of the tab-sep line 
        exclude = True
        while exclude:
            line = file.readline()
            if len(line) == 0:
                # We've reached EOF.
                self.chr = 'EOF'
                return
            line = line.rstrip()
            parts = line.split('\t')
            call = parts[column]
            callParts = call.split(':')
            #
            # The ATG field is callParts[1].  We're not interested in 0.
            exclude = callParts[1] != '1'

        # Now we have a line with a snp for this strain.
        self.chr = parts[0]
        self.position = int(parts[1])
        self.ref = parts[3]
        #The allele possibilities are comma separated. They are indexed by 
        # callParts[0].split('/')[1]
        alleles = callParts[0].split('/')
        whichAllele = int(alleles[0])-1
        try:
            self.allele = parts[4].split(',')[whichAllele]
        except IndexError:
            print >> sys.stderr, '\n\nAllele lookup failed.\n    Allele field =', parts[4], '\ncall parts =', callParts, '\nwhichAllele =', whichAllele
            sys.exit(1)

    def toResult(self):
        return Result(self.chr, self.position, self. ref, self.allele, 0,
                      self.position+1, False)

    def __str__(self):
        return self.chr + ':' + str(self.position) + ' ' + self.ref + '->' + self.allele

    def __lt__(self, indel):
        if isinstance(indel, Indel):
            if self.chr != indel.chr:
                # Currently, the VCF files SEEM TO contain chrs 1-N and X, Y, in that
                # order.  THIS PROGRAM ASSUMES THAT'S TRUE...  Hope so.
                try:
                    return int(self.chr) < int(indel.chr)
                except ValueError:
                    # One or both is an alphabetic chr name.
                    # Figure out whether either is numeric. It is "less"
                    try:
                        int(self.chr)
                        # If we get here, then we are less. (self is numeric,
                        # so indel has to be alphabetic, which we regard as
                        # greater.)
                        return True
                    except:
                        pass
                    try:
                        int(indel.chr)
                        # And if we get here, self is greater.
                        return False
                    except:
                        Pass
                    # Both were alphabetic.  Do simple string sort.
                    return self.chr < snp.chr

            # Same chr, return based on position, already ints.
            return self.position < indel.position
        print >> sys.stderr, 'Must be called with an Indel argument.'
        return NotImplemented

#Skip over vcf file headers. Return the column headers.
def skipHeaders(file):
    first = True
    while first or line[:2] == '##':
        first = False
        line = file.readline()
    if line[0:6] != '#CHROM':
        print >> sys.stderr, 'Hmmm.  Expected column headers. Exiting...'
        sys.exit(1)

    return line[1:].rstrip().split('\t')

generatedSequence = ''
def prettyizeSequence(new):
    global generatedSequence, outputFile
    flushExisting =  new == ''
    existing = generatedSequence + new
    exLen = len(existing)
    offset = 0
    while exLen >= 60:
        line = existing[offset:offset+60]
        if len(line) != 60:
            print >> sys.stderr, "ERROR: pretty has short line."
        print >> outputFile, line
        offset += 60
        exLen -= 60
    generatedSequence = existing[offset:]
    if flushExisting and len(generatedSequence) > 0:
        print >> outputFile, generatedSequence
        generatedSequence = ''

    
currentChr = ''
chrInfo = {}
offsetFile = None
outputPos = 0

def processAnEvent(ref, res, strain, dir):
    global currentChr, chrInfo, offsetFile, indelsProcessed, outputPos, outputFile

    if currentChr == '':
        offset = 0
    else:
        offset = chrInfo[currentChr]['currOffset']

    if res.chr != currentChr:
        # Set up for processing a new (maybe) chr.
        #print >> sys.stderr, 'Processing new chr', res.chr
        if offsetFile:
            # If we have an offsetFile, we have been processing a
            # chr.  Flush the previous chr's remaining sequence.
            finishChromosome(ref, currentChr, chrInfo[currentChr]['nextRefPos'])
            offsetFile.close()

        if res.chr in chrInfo:
            # We've already seen this chr!  This is a problem.
            print >> sys.stderr, 'Seeing chr', res.chr, 'again... Exiting...'
            print >> sys.stderr, 'old chr was', currentChr
            sys.exit(1)

        thisChrInfo = {}
        thisChrInfo['nextRefPos'] = 1
        thisChrInfo['currOffset'] = 0
        chrInfo[res.chr] = thisChrInfo
        currentChr = res.chr
        # Find this chr in the reference now.  This will get us the chromosome line
        # for us to print out.
        ref.findChr(currentChr)
        print >> outputFile, ref.chrLine

        outputPos = 0
        ofn = os.path.join(dir, strain + '_offsets_chr' + str(res.chr) + '.txt')
        offsetFile = open(ofn, 'w')

    # End of processing for chr changeover...

    # We have a new event, hopefully at some point past where we 
    #     previously were.
    # We need to output the ref up to this new point
    thisChrInfo = chrInfo[res.chr]

    outputError = False

    # There is at least one case of an indel being completely contained
    # in a previous del (129S1 chr1 19525839 deletes 19525840-19525879, and 
    # 19525878 deletes 19525879 )  Check for other instances of this event.
    if res.pos < int(thisChrInfo['nextRefPos']):
        # print >> sys.stderr, 'skipping event; too early: \n', res.line, '\nevent =', res.pos, 'nextRefPos =', thisChrInfo['nextRefPos']
        # Don't to anything with this event.  But before we go, check
        # for edge conditions.  If we find anything unusual, (e.g.,
        # incomplete overlap) report and exit.
        if (res.pos + (len(res.allele)-res.change)) - 1 >= res.nextRefPos:
            print >> sys.stderr, "UNUSUAL: partially pverlapping indels."
            sys.exit(2)
        return

    refSeq = ref.getRange(res.chr, thisChrInfo['nextRefPos'], res.pos-1)
    outputPos += len(refSeq)
    if len(refSeq) > 0:
        prettyizeSequence(refSeq)
    if len(res.allele) > 0:
        prettyizeSequence(res.allele)
        outputPos += len(res.allele)

    thisChrInfo['nextRefPos'] = res.nextRefPos
    if res.change != 0:
        thisChrInfo['currOffset'] += res.change
        print >> offsetFile, '%s\t%d' % (res.pos, thisChrInfo['currOffset'])

    if outputPos != (res.nextRefPos + thisChrInfo['currOffset'] - 1):
        print >> sys.stderr, '\noutputPos error: oP =', outputPos, 'res.pos:', res.pos, 'res.nextRefPos = ', res.nextRefPos, '+', thisChrInfo['currOffset']
        outputError = True

    # print >> sys.stderr, 'op:', outputPos, 'off:', thisChrInfo['currOffset']
    if outputError:
        sys.exit(1)

def processNextEvent(ref, indel, snp, strain, dir):
    processingIndel = indel < snp
    if processingIndel:
        res = indel.toResult()
    else:
        res = snp.toResult()
    processAnEvent(ref, res, strain, dir)

    return processingIndel


def processRemaining(ref, f, strainCol, strain, processingIndels, dir):
    while True:
        if processingIndels:
            event = Indel(f, strainCol)
        else:
            event = Snp(f, strainCol)
        if event.chr == 'EOF':
            # We're done.
            break
        res = event.toResult()
        processAnEvent(ref, res, strain, dir)


def finishChromosome(ref, chr, nextRefPos):
    seq = ref.getRange(chr, nextRefPos, 999999999999999)
    if len(seq) > 0:
        prettyizeSequence(seq)
    # Flush the pretty buffer.
    prettyizeSequence('')


def finishUp(ref, chrs):
    global currentChr, chrInfo, offsetFile, outputFile

    # We've finished the last events in the files.  Close out
    # the last chr we were processing.
    finishChromosome(ref, currentChr, chrInfo[currentChr]['nextRefPos'])

    # We've now output all of every chromosome that had a SNP or Indel.
    # But there could have been chromosomes without one; they aren't
    # output yet. (These may be unplaced contigs as well as complete
    # chromosomes.)
    # But they are listed in the list "chrs". So we'll output all of
    # those now.
    # FIXME DEBUG
    ref.printChrs("in Finishup")
    # The code that we're going to call in the loop below manipulated the
    # list by deleting elements.  That makes the loop unreliable, and the
    # walk misses some elements.  Make a private copy of the list, so that
    # it isn't messed with.
    chrsLocal = []
    for chr in chrs:
        chrsLocal.append(chr)
    for chr in chrsLocal:
        ref.findChr(chr)
        print >> outputFile, ref.chrLine
        seq = ref.getRange(chr, 1, 99999999999)
        if len(seq) > 0:
            prettyizeSequence(seq)
        prettyizeSequence('')

def enumerateChrs(chrs, ref, outDir):
    # It may be the case that some chromosomes don't have an indel or SNP.
    # We're outputting sequence to the new "reference" based on the processing
    # of those files.  To make sure we don't omit any chromosomes, we need
    # to know what chrs are in the reference.  So we'll make a quick scan
    # of the reference and record all the chrs we find.  To help the companion
    # program adjust_annotations.py which needs the same information, we'll 
    # also record the list in a file.
    of = open(os.path.join(outDir, 'foundChromsomes.txt'), 'w')
    for line in open(ref):
        if line[0] == '>':
            # Trim off the '>' or '>chr' and newline
            if line[:4] == '>chr':
                chr = line[4:-1]
            else:
                chr = line[1:-1]
            # Some fasta files (e.g., NCBI) have more information after the 
            # chromosome name.  Strip it off.
            chr = chr.split(' ')[0]
            chrs.append(chr)
            print >> of, chr
    of.close()

outputFile = None
def main():
    global outputFile
    opts, args = parseOptions()
    strain = args[0]

    # Track the chromosomes in this organism
    chrs = []

    # Open our files...
    ref = Reference(chrs, opts.reference)
    indels = open(opts.indels)
    snps = open(opts.snps)
    dir = opts.dir
    if opts.output:
        outputFile = open(os.path.join(dir, opts.output), 'w')
    else:
        outputFile = sys.stdout

    # First, prepare for program adjust_annotations.py by getting all the
    # chromosomes.
    enumerateChrs(chrs, opts.reference, dir)

    # Prepare for processing the indels file
    # Skip the indels vcf header cruft
    iHeaders = skipHeaders(indels)
    try:
        iStrainCol = iHeaders.index(strain)
    except ValueError:
        print >> sys.stderr, 'Hmmmm. Strain', strain, 'not found in indels header line\n', iHeaders
        sys.exit(1)
    # From here on out all lines in the indels file will be indel records.
    # Flag that we need to refresh the indel
    needIndel = True

    # Similarly set up the snps file
    sHeaders = skipHeaders(snps)
    try:
        sStrainCol = sHeaders.index(strain)
    except ValueError:
        print >> sys.stderr, 'Hmmmm. Strain', strain, 'not found in snps header line\n', sHeaders
        sys.exit(1)
    # From here on out all lines in the snps file will be snp records.
    # Flag that we need to refresh the snp
    needSnp = True
    
    # Main processing loop...
    while True:
        if needIndel:
            indel = Indel(indels, iStrainCol)
            if indel.chr == 'EOF':
                # We exhausted the indels file.  We'll clean up the snps
                # after the loop
                break
            needIndel = False
        
        if needSnp:
            snp = Snp(snps, sStrainCol)
            if snp.chr == 'EOF':
                # We exhausted the snps file.  We'll clean up the indels
                # after the loop
                break
            needSnp = False

        # Now we have an indel and a snp.  Process whichever is first.
        # This function will return True if it processed the indel,
        # False if it processed the snp.
        processedIndel = processNextEvent(ref, indel, snp, strain, dir)
        if processedIndel:
            needIndel = True
        else:
            needSnp = True

    # End of the main loop.  We have exhausted one or the other input
    # file.  Now clean up the remainder of the other file.
    if indel.chr == 'EOF':
        # Last parameter False indicates processing snps
        processRemaining(ref, snps, sStrainCol, strain, False, dir)
    elif snp.chr == 'EOF':
        processRemaining(ref, indels, iStrainCol, strain, True, dir)

    # That's about it!
    finishUp(ref, chrs)

main()
