class Reference:
    # Even though it is not coded as one (probably should fix that)
    # this is in fact a singleton.
    def __init__(self, inF):
        self.contigs = ['Y']
        self.nextChr = ''
        self.currChr = ''
        self.chrLine = ''
        self.chrSeq = ''
        self.file=open(inF)
        self.inRewindSearch = False

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
                # We'll isolate the chromosome name, and:
                #    1. If it is a contig ('NT_') stash the name for later flushing
                #    2. Treat M and MT as contig (because mm9 has it, NCBI doesn't).
                #    3. Stash the original line so we can output it.
                #
                line = line.rstrip()
                thisChr = line[1:]
                if thisChr[:3] == 'chr':
                    thisChr = thisChr[3:]
                thisChr = thisChr.split(' ')[0]
                if (thisChr[:3] == 'NT_' or thisChr[0] == 'M') and thisChr not in self.contigs :
                    self.contigs.append(thisChr)
                # Steve only wants >1
                #self.chrLine = line
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
            self.readChr(chr)
        return self.chrSeq[beg-1:end]

# End of class Reference...
