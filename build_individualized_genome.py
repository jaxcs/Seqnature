#! /usr/bin/env python
#
# build_DO_from_sanger_vcfs.py
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


# Using SNP and indel VCF files from Sanger, produce an animal-specific
# reference.
# We use a control file which describes the ranges over which a particular
# founder genome is present.
#
#   Chromosome, Position, Strain
# 
# The strain is in effect from Position to the position contained in the
# next line of the control file or to the end of the chromosome, if this
# is the last line of a chromosome (or of the file).
#

import sys
from optparse import OptionParser
from reference import Reference

# En/Dis-able some logging to stdout:
DEBUG=False
def parseOptions():
    usage = 'USAGE: %prog [options] animal-id L-R-designation\n'
    parser = OptionParser(usage=usage, version="1.0")

    parser.add_option('-c', '--control', dest='control',
                      help='path to control file')
    parser.add_option('-d', '--output-dir', dest='dir', default='.',
                      help='Directory for output files (must already '
                           'exist). (Default: current working directory.)')
    parser.add_option('-i', '--indels', dest='indels',
                      help='path to indels vcf file')
    parser.add_option('-o', '--output-sequence', dest='output',
                      help='File (fasta) into which to write the '
                           'sequence. (Default: stdout)')
    parser.add_option('-p', '--pass-only', dest='pass_only',
                      action='store_true', default=False,
                      help='Only include calls that passed the filtering '
                           'parameters used to create the vcf file.')
    parser.add_option('-r', '--reference', dest='reference',
                      help='Base name of the reference file that has been '
                           'preprocessed by reconstruct_NCBI.py, e.g., '
                           '/some/path/to/NCBIM37; not '
                           '/path/NCBIM37_ordered.fa')
    parser.add_option('-s', '--snps', dest='snps',
                      help='Path to the snps vcf file')

    (opts, args) = parser.parse_args()

    # Sanity check our options and arguments.
    # Mandatory options are reference, indels and snps
    errors = ''
    if len(args) != 2:
        errors += 'You must specify the animal ID and an L-R designation!'
    elif args[1] not in 'LR':
        errors += 'The L-R designation must be "L" or "R".'

    if not opts.control:
        errors += '\n    You must specify a control file.'
    if not opts.indels:
        errors += '\n    You must specify an indels vcf file.'
    if not opts.reference:
        errors += '\n    You must specify a reference sequence file.'
    if not opts.snps:
        errors += '\n    You must specify a snps vcf file.'

    if errors:
        parser.error(errors)
    return (opts, args)


# A common class for getting info from snps and indels.
class Result:
    def __init__(self, chr, pos, ref, allele, change,
                 nextRefPos, line, isIndel):
        self.chr = chr
        self.pos = pos
        self.ref = ref
        self.allele = allele
        self.change = change
        self.nextRefPos = nextRefPos
        self.line = line
        self.isIndel = isIndel

    def __str__(self):
        return ('I: ' if self.isIndel else 'S: ') + self.chr + ":" + \
               str(self.pos) + ' ' + self.ref + '->' + self.allele + \
               '(' + str(self.change) + ') nrp:' + str(self.nextRefPos) + \
               '\n' + self.line


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
    global opts
    strainTracker=None
    @staticmethod
    def setStrainTracker(st):
        Indel.strainTracker = st

    def __init__(self, file):
        # Get the next line in the file with a call for this strain.
        # Return a list of the tab-sep line 
        while True:
            line = file.readline()
            if len(line) == 0:
                # We've reached EOF.
                self.chr = 'EOF'
                return 
            line = line.rstrip()
            parts = line.split('\t')
            chr = parts[0]
            position = int(parts[1])
            column = Indel.strainTracker.getStrainColumn(chr, position)
            # When the strain is C57BL/6J (indicated by column number
            # -1), we don't use any indels from the Sanger file.
            # column[6] is the filter column
            if column != -1 and parts[column] != '.' and \
                    (!opts.pass or parts[6] == 'PASS')
                break

        # Now we have a line with an indel for this strain.
        self.line = line
        self.chr = chr
        ref = parts[3]
        alleleGroup = parts[4]
        call = parts[column]

        alleles = alleleGroup.split(',')
        parts1 = call.split(':')
        parts2 = parts1[0].split('/')
        if parts2[0] != parts2[1]:
            print >> sys.stderr, 'het call found at chr{0} position' \
                                 ' {1}: {2} Using first allele.' \
                                 .format(self.chr, position, call)
        #Calls use 1-based allele references.
        try:
            self.allele = alleles[int(parts2[0])-1]
        except IndexError:
            print >> sys.stderr, 'Index error processing', parts
            sys.exit(1)
            

        #
        # The way the Sanger VCFs are created, the first base in both the
        # reference and the allele will be the same. Further bases _may_
        # be the same. Make position be the first different base and fix
        # up the allele to not contain the common bases.
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
        return Result(self.chr, self.position, self.ref, self.allele,
                      self.change, self.nextRefPos, self.line, True)

    def __str__(self):
        return self.chr + ':' + str(self.position) + ' ' + self.ref + \
               '->' + self.allele + ' ' + str(self.change) + ' ' + \
               str(self.nextRefPos)

    def __lt__(self, snp):
        if isinstance(snp, Snp):
            if self.chr != snp.chr:
                # Currently, the Sanger files contain chrs 1-19 and X, in that
                # order.
                try:
                    return int(self.chr) < int(snp.chr)
                except ValueError:
                    # One is an alphabetic chr name.
                    # make sure not M or Y; these are females!
                    if self.chr[0] in 'MY':
                        print >> sys.stderr, \
                            'ERROR: Found unexpected chromosome', self.chr
                    if snp.chr[0] in 'MY':
                        print >> sys.stderr, \
                            'ERROR: Found unexpected chromosome', snp.chr
                    # One is an X  It is the greater one.
                    return snp.chr == 'X'

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
    strainTracker = None
    @staticmethod
    def setStrainTracker(st):
        Snp.strainTracker = st

    def __init__(self, file):
        global opts
        # Get the next line in the file with a call for this strain.
        # Return a list of the tab-sep line 
        while True:
            line = file.readline()
            if len(line) == 0:
                # We've reached EOF.
                self.chr = 'EOF'
                return
            line = line.rstrip()
            parts = line.split('\t')
            chr = parts[0]
            position = int(parts[1])
            column = Snp.strainTracker.getStrainColumn(chr, position)
            # When the strain is C57BL/6J (indicated by column number
            # -1), we don't use any SNPs from the Sanger file.
            # column[6] is the filter column
            if column == -1 or parts[column] == '.' or \
                    (opts.pass and parts[6] != 'PASS')
                continue
            call = parts[column]
            callParts = call.split(':')
            #
            # The ATG field is callParts[1].  We're not interested in
            # 0 or negative--only 1.  Can also be '.', in which case
            # parts is only one element long.
            if len(call_parts) > 1 and callParts[1] == '1':
                break

        # Now we have a line with a snp for this strain.
        self.line = line
        self.chr = parts[0]
        self.position = int(parts[1])
        self.ref = parts[3]
        #The allele possibilities are comma separated. They are indexed by 
        # callParts[0].split('/')[1]
        alleles = callParts[0].split('/')
        try:
            whichAllele = int(alleles[0])-1
        except:
            print >> sys.stderr, '\n\n\n\ngot allele', alleles[0], \
                'call', call, 'line', line
            sys.exit(1)
        try:
            self.allele = parts[4].split(',')[whichAllele]
        except IndexError:
            print >> sys.stderr, '\n\nAllele lookup failed.\n    ' \
                                 'Allele field =', parts[4], '\ncall ' \
                                 'parts =', callParts, '\nwhichAllele =', \
                                 whichAllele
            sys.exit(1)

    def toResult(self):
        return Result(self.chr, self.position, self. ref, self.allele, 0,
                      self.position+1, self.line, False)

    def __str__(self):
        return self.chr + ':' + str(self.position) + ' ' + self.ref + \
               '->' + self.allele

    def __lt__(self, indel):
        if isinstance(indel, Indel):
            if self.chr != indel.chr:
                # Currently, the Sanger files contain chrs 1-19 and X, in that
                # order.
                try:
                    return int(self.chr) < int(indel.chr)
                except ValueError:
                    # One is an alphabetic chr name.
                    # make sure not M or Y; these are females!
                    if self.chr[0] in 'MY':
                        print >> sys.stderr, \
                            'ERROR: Found unexpected chromosome', self.chr
                    if indel.chr[0] in 'MY':
                        print >> sys.stderr, \
                            'ERROR: Found unexpected chromosome', indel.chr
                    # One is an X  It is the greater one.
                    return indel.chr == 'X'

            # Same chr, return based on position, already ints.
            return self.position < indel.position
        print >> sys.stderr, 'Must be called with an Indel argument.'
        return NotImplemented

# End of class Snp

# Determine the strain in effect at any point
# The control file contains the starting position of each strain,
# so to know the end point of the strain, we need to read the next line.
# at the end of each chromosome and the end of the file we'll use a very
# high position to mean "until the end"
#
# Implemented this as a class because the SNP and Indel files need
# to track things separately.
class StrainTracker:
    VERY_HIGH = 99999999999
    def __init__(self, fn, header, type):
        # header is the header row from the SNP or Indels file
        StrainTracker.ordA = ord('A')
        self.type = type
        self.columns = []
        self.registerColumns(header)
        self.f = open(fn)
        self.eof = False
        # Throw away the first, header, line.
        self.f.readline()
        line = self.f.readline()
        parts = line.split(',')
        self.chr = parts[0].strip(' \t')
        self.low = int(parts[1].strip(' \t'))
        self.strain = self.strainColumn(parts[2].strip(' \t\r\n'))
        line = self.f.readline()
        parts = line.split(',')
        self.nextChr = parts[0].strip(' \t')
        if self.nextChr != self.chr:
            print >> sys.stderr, 'ERROR: not prepared for a chr change ' \
                                 'on the second line!'
            sys.exit(1)
        self.high = int(parts[1].strip(' \t'))
        self.nextStrain = self.strainColumn(parts[2].strip(' \t\r\n'))


    def nextChunk(self):
        if self.eof:
            print >> sys.stderr, 'ERROR Called nextChunk after we had ' \
                                 'already hit eof...'
            sys.exit(1)

        # Promote "next" to "current"
        if self.chr != self.nextChr:
            self.low = 1
        else:
            self.low = self.high
        self.chr = self.nextChr
        self.strain = self.nextStrain

        line = self.f.readline()
        # Special case the end of file.
        # 
        if len(line) == 0:
            self.eof = True
            self.high = StrainTracker.VERY_HIGH
            return

        # Set up "high" and "next" values
        parts = line.split(',')
        self.nextChr = parts[0].strip(' \t')
        if self.nextChr != self.chr:
            self.high = StrainTracker.VERY_HIGH
        else:
            self.high = int(parts[1].strip(' \t'))
        strainLetter = parts[2].strip(' \t\r\n')
        self.nextStrain = self.strainColumn(strainLetter)


    def registerColumns(self, header):
        #
        # The control files use the single character code for the 
        # founder strains:
        # A: A/J (Called A_J in the Sanger files)
        # B: C57BL/6J
        # C: 129S1/SvlmJ
        # D: NOD/ShiLtJ
        # E: NZO/H1LtJ
        # F: CAST/EiJ
        # G: PWK/PhJ
        # H: WSB/EiJ

        # Rather than use a hash, we're going to index into a table
        # using the ordinals of the code letters.
        strainNames = ['A_J', 'C57BL', '129S1', 'NOD', 'NZO', 'CAST',
                       'PWK', 'WSB']
        for n in range(len(strainNames)):
            self.columns.append(header.index(strainNames[n]))
        self.columns[1] = -1 # FLAG: for C57BL/6J, we're using the reference,
                        # not the Sanger C57BL, which is C57BL/6NJ

    def strainColumn(self, strain):
        try:
            col =  self.columns[ord(strain) - StrainTracker.ordA]
        except:
            print >> sys.stderr, 'range exception', strain, ord(strain), \
                ord(strain) - StrainTracker.ordA, self.columns
            sys.exit(1)
        return col
    
    def getStrainColumn(self, chr, pos):
        while chr != self.chr:
            self.nextChunk()
        if pos < self.low:
            # I think that this means something is out of order.
            # Report it and bail.
            print >> sys.stderr, "Pos below current low point. " \
                                 "Ordering problem?"
            sys.exit(1)
        while pos >= self.high:
            self.nextChunk()
        return self.strain

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

#
# Sometimes an animal will have a chromosome that is entirely B6.
# This means that there would have been no events in that chromosome,
# so processAnEvent() would merrily go from chr10 to, say, chr12.
# This routine helps us catch that event, and outputs the skipped
# chromosome(s).
#

# Include '' to handle the initial case. Include a tombstone/known
# elephant 'Done' to handle emitting X.
chrsInOrder=['', '1', '2', '3', '4', '5', '6', '7', '8','9', '10', '11',
             '12', '13', '14', '15', '16', '17', '18', '19', 'X', 'Done']

def skippedChrs(old, new, ref, lr):
    oIdx = chrsInOrder.index(old)
    nIdx = chrsInOrder.index(new)
    if nIdx == oIdx + 1:
        return
    #
    # So we've skipped at least one chromosome.
    # Output it/them.
    for n in range(oIdx+1, nIdx):
        # output the segment header, followed by the entire chr.
        print ">" + chrsInOrder[n] + lr
        emitEntireChromosome(ref, chrsInOrder[n])

    # There... Doesn't that feel better?


def processAnEvent(ref, res, animalId):
    global currentChr, chrInfo, offsetFile, indelsProcessed, outputPos, \
        lr, DEBUG

    if currentChr == '':
        offset = 0
    else:
        offset = chrInfo[currentChr]['currOffset']

    if DEBUG:
        print >> sys.stderr, 'op:', outputPos, 'off:', offset, res,

    if res.chr != currentChr:
        # Set up for processing a new (maybe) chr.
        #print >> sys.stderr, 'Processing new chr', res.chr
        if offsetFile:
            # If we have an offsetFile, we have been processing a
            # chr.  Flush the previous chr's remaining sequence.
            finishChromosome(ref, currentChr,
                             chrInfo[currentChr]['nextRefPos'])
            offsetFile.close()

        if res.chr in chrInfo:
            # We've already seen this chr!  This is a problem.
            print >> sys.stderr, 'Seeing chr', res.chr, 'again... Exiting...'
            print >> sys.stderr, 'old chr was', currentChr
            sys.exit(1)

        # Check whether we've skipped a chromosome or two...
        skippedChrs(currentChr, res.chr, ref, lr)

        mode = 'w'
        thisChrInfo = {}
        thisChrInfo['nextRefPos'] = 1
        thisChrInfo['currOffset'] = 0
        chrInfo[res.chr] = thisChrInfo
        currentChr = res.chr
        # Find this chr in the reference now.  This will get us the
        # chromosome line for us to print out.
        ref.findChr(currentChr)
        print ref.chrLine + lr

        outputPos = 0
        offsetFile = open(animalId + '_offsets_chr' + str(res.chr) +
                          '.txt', mode)

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
        print >> sys.stderr, 'skipping event; too early: \n', res.line, \
            '\nevent =', res.pos, 'nextRefPos =', thisChrInfo['nextRefPos']
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
        print >> sys.stderr, '\noutputPos error: oP =', outputPos, \
            'res.pos:', res.pos, 'res.nextRefPos = ', res.nextRefPos, '+', \
            thisChrInfo['currOffset']
        outputError = True

    #print >> sys.stderr, 'op:', outputPos, 'off:', thisChrInfo['currOffset']
    if outputError:
        sys.exit(1)

def processNextEvent(ref, indel, snp, animalId):
    processingIndel = indel < snp
    if processingIndel:
        res = indel.toResult()
    else:
        res = snp.toResult()
    processAnEvent(ref, res, animalId)

    return processingIndel


def processRemaining(ref, f, animalId, processingIndels):
    while True:
        if processingIndels:
            event = Indel(f)
        else:
            event = Snp(f)
        if event.chr == 'EOF':
            # We're done.
            break
        res = event.toResult()
        processAnEvent(ref, res, animalId)

def emitEntireChromosome(ref, chr):
    seq = ref.getRange(chr, 1, 999999999999999)
    prettyizeSequence(seq)
    prettyizeSequence('')

def finishChromosome(ref, chr, nextRefPos):
    seq = ref.getRange(chr, nextRefPos, 999999999999999)
    if len(seq) > 0:
        prettyizeSequence(seq)
    # Flush the pretty buffer.
    prettyizeSequence('')


def finishUp(ref, reference, male, lr):
    global currentChr, chrInfo, offsetFile, outputFile

    # We've finished the last events in the files.  Close out
    # the last chr we were processing.
    finishChromosome(ref, currentChr, chrInfo[currentChr]['nextRefPos'])

    # Make sure that we finished all the chromosomes...
    skippedChrs(currentChr, 'Done', ref, lr)

    # We're done with the contigs for which Sanger had variant information.
    # If this is the 'R' pass, we need to append the unplaced contigs ,
    #("NT_*") MT, and Y if this sample is male.
    if lr == 'L':
        return
    if male:
        for line in open(reference + '_Y.fa'):
            print >> outputFile, line,
    for line in open(reference + '_contigs.fa'):
        print >> outputFile, line,

lr = ''
outputFile = sys.stdout
opts = None
def main():
    global lr, outputFile, opts

    opts, args = parseOptions()
    lr = args[1]
    animalId = args[0] + '_' + lr
    male = animalId[0] == 'M'

    # Open our files...
    ref = Reference(opts.reference + '_ordered.fa')
    indels = open(opts.indels)
    snps = open(opts.snps)

    dir = opts.dir
    
    if opts.output:
        outputFile = open(os.path.join(dir, opts.output), 'w')
    else:
        outputFile = sys.stdout
    

    # Prepare for processing the indels file
    # Skip the indels vcf header cruft
    iHeaders = skipHeaders(indels)
    iStrainTracker = StrainTracker(opts.control, iHeaders, 'I')
    Indel.setStrainTracker(iStrainTracker)

    # From here on out all lines in the indels file will be indel records.
    # Flag that we need to refresh the indel
    needIndel = True

    # Similarly set up the snps file
    sHeaders = skipHeaders(snps)
    sStrainTracker = StrainTracker(opts.control, sHeaders, 'S')
    Snp.setStrainTracker(sStrainTracker)

    # From here on out all lines in the snps file will be snp records.
    # Flag that we need to refresh the snp
    needSnp = True
    
    # Main processing loop...
    while True:
        if needIndel:
            indel = Indel(indels)
            if indel.chr == 'EOF':
                # We exhausted the indels file.  We'll clean up the snps
                # after the loop
                break
            needIndel = False
        
        if needSnp:
            snp = Snp(snps)
            if snp.chr == 'EOF':
                # We exhausted the snps file.  We'll clean up the indels
                # after the loop
                break
            needSnp = False

        # Now we have an indel and a snp.  Process whichever is first.
        # This function will return True if it processed the indel,
        # False if it processed the snp.
        processedIndel = processNextEvent(ref, indel, snp, animalId)
        if processedIndel:
            needIndel = True
        else:
            needSnp = True

    # End of the main loop.  We have exhausted one or the other input
    # file.  Now clean up the remainder of the other file.
    if indel.chr == 'EOF':
        # Last parameter False indicates processing snps
        processRemaining(ref, snps, animalId, False)
    elif snp.chr == 'EOF':
        processRemaining(ref, indels, animalId, True)

    # That's about it!
    finishUp(ref, opts.reference, male, lr)

main()
