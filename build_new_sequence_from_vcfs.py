#! /usr/bin/env python
#
# build_reference_from_sanger_vcfs.py
#

# Using SNP and indel VCF files from Sanger, produce pseudo-reference
# genomes for any strain represented in the VCF.
#

#   Author: Al Simons 
#   Copyright (C) 2012, 2015  The Jackson Laboratory
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
#   You should have received a copy of the GNU Affero General Public
#   License along with this program.  If not, see
#   <http://www.gnu.org/licenses/>.

import sys
import os
from optparse import OptionParser
from pyfaidx import Fasta


def parse_options():
    usage = 'USAGE: %prog [options] strain\n    Strain must be a ' \
            'column header in the two VCF files.'

    parser = OptionParser(usage=usage, version="1.2")

    parser.add_option('-d', '--output-dir', dest='dir', default='.',
                      help='Directory for output files (must already '
                           'exist). (Default: current working '
                           'directory.)')
    parser.add_option('-i', '--indels', dest='indels',
                      help='path to indels vcf file')
    parser.add_option('-o', '--output',
                      help='File into which the genome is written. '
                           '(Default: stdout)')
    parser.add_option('-p', '--pass-only',
                      action='store_true', default=False,
                      help='Only include calls marked as PASS in the '
                           'FILTER column of the vcf file.')
    parser.add_option('-r', '--reference',
                      help='Path to the reference file used as the '
                           'base.')
    parser.add_option('-s', '--snps',
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
    return opts, args

het_file = None


def log_het_call(msg):
    global het_file, args, opts
    if het_file is None:
        fn = os.path.join(opts.dir, args[0] + '_het_calls.log')
        het_file = open(fn, 'w')
    print >> het_file, msg


class Reference:
    # Even though it is not coded as one this is in fact a singleton.

    def __init__(self, chrs, in_f):
        # pyfaidx uses (as the name would suggest), a .fa.fai index
        # file.  If not present in the same directory as the .fa
        # file, it will be generated the first time the .fa is opened.
        # Therefore, on first use, the next line will take a while.
        self.fasta = Fasta(in_f)
        self.chrs = chrs

    def get_range(self, chr, beg, end):
        return self.fasta[chr][beg - 1:end]

# End of class Reference...


# A common class for getting info from snps and indels.
class Result:
    def __init__(self, chr, pos, ref, allele, change, next_ref_pos,
                 is_indel):
        self.chr = chr
        self.pos = pos
        self.ref = ref
        self.allele = allele
        self.change = change
        self.next_ref_pos = next_ref_pos
        self.is_indel = is_indel

    def __str__(self):
        return ('I: ' if self.is_indel else 'S: ') + self.chr + ":" + \
            str(self.pos) + ' ' + self.ref + '->' + self.allele + \
            '(' + str(self.change) + ') nrp:' + \
            str(self.next_ref_pos)


#
# A class to track indels.  Use toResult() to get results.
# Members:
#     chr - self explanatory
#     position - This is NOT the position contained in the VCF file,
#                which was the last common base between the reference
#                and the allele.  This is some number of bases later,
#                of the first different base.  Consider the case of
#                the position chr1 3000742, ref: GT, alleles
#                G,GTT,GGTTTT.
#                The first allele is a one base deletion.  The second is
#                a one base insertion, 2 bases after the listed position
#                the third is a change of the second base and then an
#                insertion of 4 bases. Position for this allele would
#                be one after the pos in the file.
#     allele -   The bases (for an insert) starting at the first
#                different base. Null string for a deletion.
#     change -   the effect of this allele on the new coordinates.
#                len(allele) - len(ref)
#     nextRefPos - The next base position that is to be taken from the
#                reference.
#
class Indel:
    def __init__(self, file, haplotype_column, i_columns):
        global opts
        # Get the next line in the file with a call for this strain.
        # Return a list of the tab-sep line 
        # Column 6 (the 7th column) is the FILTER column, which records
        # whether or not the call passed all the filtering criteria
        # at the time the VCF file was created.
        first = True
        short = False
        parts = []
        ref_call = False
        while first or short or \
                (opts.pass_only and
                    parts[i_columns['FILTER']] != 'PASS') or \
                parts[haplotype_column] == '.' or \
                ref_call:
            line = file.readline()
            if len(line) == 0:
                # We've reached EOF.
                self.chr = 'EOF'
                return
            line = line.rstrip()
            parts = line.split('\t')
            short = len(parts) < 9
            first = False

            # Find out if this is a ref call for this strain.
            # Since in the event of a het call, we arbitrarily use the
            # first called allele, we'll determine that this is a ref
            # call if the first identified allele is 0.  However, if it
            # is het, we'll report it.

            self.chr = parts[i_columns['#CHROM']]
            position = int(parts[i_columns['POS']])
            call = parts[haplotype_column]

            parts1 = call.split(':')
            parts2 = parts1[0].split('/')

            # Is this a ref call? Continue if so, by setting up the
            # loop condition
            ref_call = parts2[0] == '0'

            # Report a het call.
            if parts2[0] == '0' and parts2[0] != parts2[1]:
                log_het_call('Indel het REF call at chr{0}'
                    'position {1}. Using ref; continuing.'.format(
                        self.chr, position))

        # Now we have a line with an indel for this strain.
        ref = parts[i_columns['REF']]
        allele_group = parts[i_columns['ALT']]
        alleles = allele_group.split(',')
        #Calls use 1-based allele references.
        try:
            self.allele = alleles[int(parts2[0]) - 1]
        except IndexError:
            print >> sys.stderr, 'Index error processing', parts
            sys.exit(1)

        # Warn about a het call
        if parts2[0] != parts2[1]:
            log_het_call('Indel het call found at chr{0} position '
                '{1}: {2} Using alt allele {3} ({4}).'.format(
                    self.chr, position, call, parts2[0],
                    self.allele))
        #
        # The way the Sanger VCFs are created, the first base in both
        # the reference and the allele will be the same. Further bases
        # _may_ be the same. Make position be the first different base
        # and fix up the allele to not contain the common bases.
        allele = self.allele
        min_len = min(len(ref), len(allele))
        diff_pos = 0
        for n in range(min_len):
            if ref[n] != allele[n]:
                diff_pos = n
                break

        # If diff_pos == 0, the reference and allele matched throughout
        # their length.  This is either an insertion after the
        # reference, or a deletion.
        # Either way, the difference position is the base after the
        # shorter of the two.
        if diff_pos == 0:
            diff_pos = min_len
        self.position = position + diff_pos
        self.allele = self.allele[diff_pos:]
        self.ref = ref[diff_pos:]

        # A somewhat common pattern is for a few bases at the beginning
        # to be deleted, followed by many bases that are the same.
        # This is necessitated by the multi-strain nature of the file. 
        # Some other strain may have deleted many more bases here.
        len_diff = len(self.ref) - len(self.allele)
        if len_diff > 0 and self.ref[len_diff:] == self.allele:
            self.ref = self.ref[:len_diff]
            self.allele = ''

        # Also, the change in position will be len(allele) - len(ref).
        # This will be negative for a deletion, and positive for an
        # insertion.
        self.change = len(self.allele) - len(self.ref)
        self.next_ref_pos = self.position + len(self.ref)

    def to_result(self):
        return Result(self.chr, self.position, self.ref, self.allele,
                      self.change, self.next_ref_pos, True)

    def __str__(self):
        return self.chr + ':' + str(self.position) + ' ' + self.ref + \
            '->' + self.allele + ' ' + str(self.change) + ' ' + \
            str(self.next_ref_pos)

    def __lt__(self, snp):
        if isinstance(snp, Snp):
            if self.chr != snp.chr:
                # Currently, the VCF files SEEM TO contain chrs 1-N
                # and X, Y, in that order.
                try:
                    return int(self.chr) < int(snp.chr)
                except ValueError:
                    # One or both is an alphabetic chr name.
                    # Figure out whether either is numeric. It is "less"
                    try:
                        int(self.chr)
                        # If we get here, then we are less. (self is
                        # numeric, so snp has to be alphabetic, which
                        # we regard as greater.)
                        return True
                    except:
                        pass
                    try:
                        int(snp.chr)
                        # And if we get here, self is greater.
                        return False
                    except:
                        pass
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
    def __init__(self, file, haplotype_column, s_columns):
        global opts
        # Get the next line in the file with a call for this strain.
        # Return a list of the tab-sep line 
        # Column 6 is the FILTER column.  See Indel.__init__().
        exclude = True
        parts = []
        call_parts = []
        ref_call = False
        while exclude or \
              (opts.pass_only and
                  parts[s_columns['FILTER']] != 'PASS') or \
              ref_call:
            line = file.readline()
            if len(line) == 0:
                # We've reached EOF.
                self.chr = 'EOF'
                return
            line = line.rstrip()
            parts = line.split('\t')
            if len(parts) < 9:
                continue
            call = parts[haplotype_column]
            call_parts = call.split(':')
            # Old style Sanger VCFs indicated no call for this
            # haplotype by a dot. Newer forms indicate by calling the
            # reference with '0/0'.
            exclude = len(call_parts) < 2

            self.chr = parts[0]
            self.position = int(parts[1])
            alleles = call_parts[0].split('/')

            # Is this a ref call? Continue if so, by setting up the
            # loop condition
            ref_call = alleles[0] == '0'

            # Report a het call.
            if alleles[0] == '0' and alleles[0] != alleles[1]:
                log_het_call('SNP het REF call at chr{0}'
                    'position {1}. Using ref; continuing.'.format(
                        self.chr, self.position))

        # Now we have a line with a snp for this strain.
        #print >> sys.stderr, "SNP-Found", call, line
        self.ref = parts[3]
        # The allele possibilities are comma separated. They are indexed
        # by call_parts[0].split('/')[1]
        #
        which_allele = int(alleles[0]) - 1

        try:
            self.allele = parts[4].split(',')[which_allele]
        except IndexError:
            print >> sys.stderr, '\n\nAllele lookup failed.\n    ' \
                                 'Allele field =', parts[4], \
                                 '\ncall parts =', call_parts, \
                                 '\nwhich allele =', which_allele
            sys.exit(1)

        # Warn about a het call
        if alleles[0] != alleles[1]:
            log_het_call('SNP het call found at chr{0} position '
                '{1}: {2} Using alt allele {3} ({4}).'.format(
                    self.chr, self.position, call, which_allele,
                    self.allele))

    def to_result(self):
        return Result(self.chr, self.position, self. ref, self.allele,
                    0, self.position + 1, False)

    def __str__(self):
        return self.chr + ':' + str(self.position) + ' ' + self.ref + \
            '->' + self.allele

    def __lt__(self, indel):
        if isinstance(indel, Indel):
            if self.chr != indel.chr:
                # Currently, the VCF files SEEM TO contain chrs 1-N
                # and X, Y, in that order.
                try:
                    return int(self.chr) < int(indel.chr)
                except ValueError:
                    # One or both is an alphabetic chr name.
                    # Figure out whether either is numeric. It is "less"
                    try:
                        int(self.chr)
                        # If we get here, then we are less. (self is
                        # numeric, so indel has to be alphabetic, which
                        # we regard as greater.)
                        return True
                    except:
                        pass
                    try:
                        int(indel.chr)
                        # And if we get here, self is greater.
                        return False
                    except:
                        pass
                    # Both were alphabetic.  Do simple string sort.
                    return self.chr < indel.chr

            # Same chr, return based on position, already ints.
            return self.position < indel.position
        print >> sys.stderr, 'Must be called with an Indel argument.'
        return NotImplemented


#Skip over vcf file headers. Return the column headers.
def skip_headers(f):
    line = '##'
    while line[:2] == '##':
        line = f.readline()
    if line[:6] != '#CHROM':
        print >> sys.stderr, 'Hmmm.  Expected column headers; ' \
                             'got {0}.\nExiting...'.format(line)
        sys.exit(1)
    columns = {}
    parts = line.split('\t')
    for n in range(len(parts)):
        columns[parts[n].strip()] = n

    return columns, line

generated_sequence = ''


def prettyize_sequence(new):
    global generated_sequence, output_file
    flush_existing = new == ''
    existing = generated_sequence + str(new)
    ex_len = len(existing)
    offset = 0
    while ex_len >= 60:
        line = existing[offset:offset + 60]
        if len(line) != 60:
            print >> sys.stderr, "ERROR: pretty has short line."
        print >> output_file, line
        offset += 60
        ex_len -= 60
    generated_sequence = existing[offset:]
    if flush_existing and len(generated_sequence) > 0:
        print >> output_file, generated_sequence
        generated_sequence = ''


current_chr = ''
chr_info = {}
offset_file = None
output_pos = 0


def process_an_event(ref, res, strain, dir):
    global current_chr, chr_info, offset_file, \
        output_pos, output_file

    if res.chr != current_chr:
        # Set up for processing a new (maybe) chr.
        #print >> sys.stderr, 'Processing new chr', res.chr
        if offset_file:
            # If we have an offsetFile, we have been processing a
            # chr.  Flush the previous chr's remaining sequence.
            finish_chromosome(ref, current_chr, chr_info[current_chr]
                             ['nextRefPos'])
            offset_file.close()

        if res.chr in chr_info:
            # We've already seen this chr!  This is a problem.
            # We expect the indel and snp files to have contiguous
            # chromosomes.
            print >> sys.stderr, 'Seeing chr {0} again... Exiting...'\
                '\n    old chr was {1}'.format(res.chr, current_chr)
            sys.exit(1)

        this_chr_info = {}
        this_chr_info['nextRefPos'] = 1
        this_chr_info['currOffset'] = 0
        chr_info[res.chr] = this_chr_info
        current_chr = res.chr
        # Find this chr in the reference now.  This will get us the
        # chromosome line for us to print out.

        ln = ref.fasta[current_chr].long_name

        # FIXME  Work around a bug in pyfaidx, seen 20150401:
        # If the last line of a chr is a full 60 characters,
        # the long name of the following chr will be missing
        # its first character.  There are two instances of
        # this in GRCm38_68.fa, of which chr 3, affecting
        # chr 4, is of interest to us (the other is a scaffold,
        # which we aren't processing.  The bug fix, below,
        # will catch these cases, and will continue to be OK
        # when the bug is fixed.  Hopefully, we can pull it
        # out at some point.
        if current_chr[0] != ln[0]:
            ln = current_chr[0] + ln
        # FIXME: End of bug workaround.

        print >> output_file, '>{0}'.format(ln)

        output_pos = 0
        ofn = os.path.join(dir, '{0}_offsets_chr{1}.txt'.format(
            strain, current_chr))
        offset_file = open(ofn, 'w')

    # End of processing for chr changeover...

    # We have a new event, hopefully at some point past where we 
    #     previously were.
    # We need to output the ref up to this new point
    this_chr_info = chr_info[res.chr]

    output_error = False

    # There is at least one case of an indel being completely contained
    # in a previous del (129S1 chr1 19525839 deletes 19525840-19525879,
    # and 19525878 deletes 19525879 )  Check for other instances of
    # this event.
    if res.pos < int(this_chr_info['nextRefPos']):
        # print >> sys.stderr, 'skipping event; too early: \n', \
        # res.line, '\nevent =', res.pos, 'nextRefPos =', \
        # this_chr_info['nextRefPos']
        # Don't to anything with this event.  But before we go, check
        # for edge conditions.  If we find anything unusual, (e.g.,
        # incomplete overlap) report and exit.
        if (res.pos + (len(res.allele) - res.change)) - 1 >= \
                res.next_ref_pos:
            print >> sys.stderr, \
                "UNUSUAL: partially overlapping indels."
            sys.exit(2)
        return

    ref_seq = ref.get_range(res.chr, this_chr_info['nextRefPos'],
                            res.pos - 1)
    output_pos += len(ref_seq)
    if len(ref_seq) > 0:
        prettyize_sequence(ref_seq)
    if len(res.allele) > 0:
        prettyize_sequence(res.allele)
        output_pos += len(res.allele)

    this_chr_info['nextRefPos'] = res.next_ref_pos
    if res.change != 0:
        this_chr_info['currOffset'] += res.change
        print >> offset_file, '%s\t%d' % (res.pos,
                                         this_chr_info['currOffset'])

    if output_pos != (res.next_ref_pos +
                    this_chr_info['currOffset'] - 1):
        print >> sys.stderr, '\noutputPos error: oP =', output_pos, \
            'res.pos:', res.pos, 'res.nextRefPos = ', \
            res.next_ref_pos, '+', this_chr_info['currOffset']
        output_error = True

    # print >> sys.stderr, 'op:', outputPos, 'off:', \
    #     this_chr_info['currOffset']
    if output_error:
        sys.exit(1)


def process_next_event(ref, indel, snp, strain, dir):
    processing_indel = indel < snp
    if processing_indel:
        res = indel.to_result()
    else:
        res = snp.to_result()
    process_an_event(ref, res, strain, dir)

    return processing_indel


def process_remaining(ref, f, strain_col, strain,
                      processing_indels, dir, columns):
    if processing_indels:
        print 'Processing remaining indels'
    else:
        print 'Processing remaining snps'
    while True:
        if processing_indels:
            event = Indel(f, strain_col, columns)
        else:
            event = Snp(f, strain_col, columns)
        if event.chr == 'EOF':
            # We're done.
            break
        res = event.to_result()
        process_an_event(ref, res, strain, dir)


def finish_chromosome(ref, chr, next_ref_pos):
    print 'Finishing chr', chr
    seq = ref.get_range(chr, next_ref_pos, sys.maxint)
    if len(seq) > 0:
        prettyize_sequence(seq)
    # Flush the pretty buffer.
    prettyize_sequence('')


def finish_up(ref, chrs):
    global current_chr, chr_info, offset_file, output_file

    # We've finished the last events in the files.  Close out
    # the last chr we were processing.
    finish_chromosome(ref, current_chr,
                      chr_info[current_chr]['nextRefPos'])

    # We've now output all of every chromosome that had a SNP or Indel.
    # But there could have been chromosomes without one; they aren't
    # output yet. (These may be unplaced contigs as well as complete
    # chromosomes.)
    # But they are listed in the list "chrs". So we'll output all of
    # those now.

    # The global chr_info contains a key for each chr we've processed.
    # Make a new list containing only those chrs we haven't seen, and
    # write them out unchanged.

    chrs_local = []
    for chr in chrs:
        if chr not in chr_info:
            chrs_local.append(chr)

    for chr in chrs_local:
        print 'Writing out', chr

        print >> output_file, '>' + ref.fasta[chr].long_name
        seq = ref.get_range(chr, 1, sys.maxint)
        if len(seq) > 0:
            prettyize_sequence(seq)
        prettyize_sequence('')


def enumerate_chrs(chrs, ref, out_dir):
    # It may be the case that some chromosomes don't have an indel or
    # SNP. We're outputting sequence to the new "reference" based on
    # the processing of those files.  To make sure we don't omit any
    # chromosomes, we need to know what chrs are in the reference.
    # Fortunately, pyfaidx gives us an easy way to do that.
    # To help the companion program adjust_annotations.py
    # which needs the same information, we'll
    # also record the list in a file.

    # This is also where we drop out all the scaffolds included by the
    # GRC, if any.  We don't want to process those.
    of = open(os.path.join(out_dir, 'foundChromsomes.txt'), 'w')
    c = ref.fasta.keys()
    for chr in sorted(c):
        # Filter out unplaced extents / scaffolds.
        if len(chr) < 4:
            chrs.append(chr)
            print >> of, chr
    of.close()

output_file = None
opts = None
args = None


def main():
    global output_file, opts, args
    opts, args = parse_options()
    strain = args[0]

    # Track the chromosomes in this organism
    chrs = []

    # Open our files...
    ref = Reference(chrs, opts.reference)
    indels = open(opts.indels)
    snps = open(opts.snps)
    dir = opts.dir
    indel = None
    snp = None

    if opts.output:
        output_file = open(os.path.join(dir, opts.output), 'w')
    else:
        output_file = sys.stdout

    # First, prepare for program adjust_annotations.py by getting all
    # the chromosomes.
    enumerate_chrs(chrs, ref, dir)

    # Prepare for processing the indels file
    # Skip the indels vcf header cruft
    i_columns, i_header_line = skip_headers(indels)
    try:
        i_strain_col = i_columns[strain]
    except KeyError:
        print >> sys.stderr, 'Hmmmm. Strain', strain, \
            'not found in indels header line\n', i_header_line
        sys.exit(1)
    # From here on out all lines in the indels file will be indel
    # records. Flag that we need to refresh the indel
    need_indel = True

    # Similarly set up the snps file
    s_columns, s_header_line = skip_headers(snps)
    try:
        s_strain_col = s_columns[strain]
    except KeyError:
        print >> sys.stderr, 'Hmmmm. Strain', strain, \
            'not found in snps header line\n', s_header_line
        sys.exit(1)
    # From here on out all lines in the snps file will be snp records.
    # Flag that we need to refresh the snp
    need_snp = True

    # Main processing loop...
    while True:
        if need_indel:
            indel = Indel(indels, i_strain_col, i_columns)
            if indel.chr == 'EOF':
                # We exhausted the indels file.  We'll clean up the snps
                # after the loop
                break
            need_indel = False

        if need_snp:
            snp = Snp(snps, s_strain_col, s_columns)
            if snp.chr == 'EOF':
                # We exhausted the snps file.  We'll clean up the indels
                # after the loop
                break
            need_snp = False

        # Now we have an indel and a snp.  Process whichever is first.
        # This function will return True if it processed the indel,
        # False if it processed the snp.
        processed_indel = process_next_event(ref, indel, snp,
                                             strain, dir)
        if processed_indel:
            need_indel = True
        else:
            need_snp = True

    # End of the main loop.  We have exhausted one or the other input
    # file.  Now clean up the remainder of the other file.
    if indel.chr == 'EOF':
        # Last parameter False indicates processing snps
        process_remaining(ref, snps, s_strain_col, strain, False, dir,
                          s_columns)
    elif snp.chr == 'EOF':
        process_remaining(ref, indels, i_strain_col, strain, True, dir,
                          i_columns)

    # That's about it!
    finish_up(ref, chrs)

if __name__ == '__main__':
    main()
