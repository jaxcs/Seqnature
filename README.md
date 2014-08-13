#Seqnature

If you use Seqnature for your published research, please cite:

Munger SC, Raghupathy N, Choi KB, Simons AK, Gatti DM, Hinerfeld DA, Svenson KL, 
Keller MP, Attie AD, Hibbs MA, Graber JH, Chesler EJ, and Churchill GA. 
RNA-seq alignment to individualized genomes improves transcript abundance 
estimates in multiparent populations. GENETICS, In Press. 

The programs' author is Al Simons (al.simons at jax.org).  Please contact him
to report any problems with the programs.

_Seqnature_ is a set of programs to create custom sequences from a reference 
sequence and a pair of VCF files containing SNP and Indel information.  
It was created for, and most heavily tested on, mouse.

If you have gtf files with annotations based on the same reference, the
programs can create an updated gtf, reflecting the Indels.

##Creating new sequence and annotations for strains
This section discusses the programs used to create strain-specific references
and annotation files, __build_new_sequence_from_vcfs.py__ and
__adjust_annotations.py__.  For information about building individualized
reference sequences and annotation files, please see the section
"Creating new sequence and annotations for individuals", below.

The program __build_new_sequence_from_vcfs.py__ takes a reference sequence, 
and two vcf files, one each of SNPs and Indels.  It produces a new
sequence file and a set of bookkeeping files that are used by the second
program in the set, __adjust_annotations.py__.  The following options are
required:
* --indels
* --reference
* --snps

The program __adjust_annotations.py__ takes a gtf file with coordinates matching
the reference sequence previously processed, and the bookkeeping files produced
by __build_new_sequence_from_vcfs.py__.  It produces an updated gtf file reflecting
the effects of the Indels.  While most commonly used on gtf files, it is able 
to process any tab separated file with columns containing chromosome and 
positions within the chromosome that are to be updated.

The _Seqnature_ programs expect chromosome names of the form "1", "2", "X", etc,
not "chr1", "chr2", etc.

###Using build_new_sequence_from_vcfs.py

    $ ./build_new_sequence_from_vcfs.py --help
    Usage: build_new_sequence_from_vcfs.py [options] strain
        Strain must be a column header in the two VCF files.

    Options:
      --version             show program's version number and exit
      -h, --help            show this help message and exit
      -d DIR, --output-dir=DIR
                            Directory for output files (must already exist).
                            (Default: current working directory.)
      -i INDELS, --indels=INDELS
                            path to indels vcf file
      -o OUTPUT, --output=OUTPUT
                            File into which the genome is written. (Default:
                            stdout)
      -p, --pass-only       Only use SNPs and Indels that passed quality
                            filters
      -r REFERENCE, --reference=REFERENCE
                            Path to the reference file used as the base.
      -s SNPS, --snps=SNPS  Path to the snps vcf file

###Using adjust_annotations.py

    $ ./adjust_annotations.py --help
    Usage: adjust_annotations.py [options] annotation-file strain-identifier

    Options:
      --version             show program's version number and exit
      -h, --help            show this help message and exit
      -c CHRCOL, --chr-column=CHRCOL
                            Column with chromosome info. (Required, one-based.)
      -C COMMENTS, --comments-file=COMMENTS
                            File to write adjustment comments to. (Default: no
                            comments file.)
      -d DIR, --directory=DIR
                            path to directory containing offset files (default:
                            current directory)
      -e ENDCOL, --end-column=ENDCOL
                            The column containing the position of the end of the
                            feature. (Required, one-based.)
      -n NEAR, --near=NEAR  Report if a feature position is within NEAR bases of
                            an indel. (Optional, default is no report.)
      -o OFN, --output-file=OFN
                            Output filename (Optional, default: stdout)
      -p POSCOLS, --position-columns=POSCOLS
                            Comma-separated list of other columns with positions
                            to be adjusted. (Optional, one-based.)
      -q, --quiet           Operate quietly. Do not report chromosomes and contigs
                            not in the Indel vcf file.
      -s STARTCOL, --start-column=STARTCOL
                            The column containing the position of the end of the
                            feature. (Required, one-based.)
      -t TYPECOL, --type-column=TYPECOL
                            Column with feature type info. (Optional, one-based.)

##Creating new sequence and annotations for individuals
This section describes the programs used to create new reference sequences and 
annotation files for individual animals which have been haplotyped: 
__reconstruct_NCBI.py__, __build_individualized_genome.py__ and
__adjust_individualized_annotations.py__.

An example shell script for running all of these steps other than reordering the
reference, is included as do_DO.sh.

###Reordering the reference for performance
The program __reconstruct_NCBI.py__ reorders the original genome reference into 
chromosome order, and splits out all the NT_ unplaced extents and chromosomes Y 
and MT.  This allows marked improvements in the performance of the program that 
builds the individualized genome.  This program only has to be run once, not
once per animal.

###Building the individualized genome reference
After the reference has been reordered, the program __build_individualized_genome.py__
will create one chromatid's reference for the individual.  This program is run twice,
once per haplotype file ('L' and 'R').  At the end of the 'R' run, the program
includes the unplaced NT_ contigs
and chr MT from the original reference.  If the sample is male, 
the unmodified chromosome Y is appended as well.  After running the program twice,
the output files are concatenated to form the complete diploid reference.

###Adjusting the annotations file
Once the genome has been created by two invocations of __build_individualized_genome.py__,
the annotation gtf file can be adjusted by program __adjust_individualized_annotations.py__.
Just as with __build_individualized_genome.py__, __adjust_individualized_annotations.py__ is
run twice, once per haplotype block, and the output from the two runs are concatenated.

###Using reconstruct_NCBI.py
Running the program without any arguments produces the following usage information:

    $ ./reconstruct_NCBI.py
    
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

###Using build_individualized_genome.py
    ./build_individualized_genome.py --help
    Usage: build_individualized_genome.py [options] animal-id L-R-designation

    Options:
      --version             show program's version number and exit
      -h, --help            show this help message and exit
      -c CONTROL, --control=CONTROL
                            path to control file
      -d DIR, --output-dir=DIR
                            Directory for output files (must already exist).
                            (Default: current working directory.)
      -i INDELS, --indels=INDELS
                            path to indels vcf file
      -o OUTPUT, --output-sequence=OUTPUT
                            File (fasta) into which to write the sequence.
                            (Default: stdout)
      -p, --pass-only       Only use SNPs and Indels that passed quality
                            filters
      -r REFERENCE, --reference=REFERENCE
                            Base name of the reference file that has been
                            preprocessed by reconstruct_NCBI.py, e.g.,
                            /some/path/to/NCBIM37; not /path/NCBIM37_ordered.fa
      -s SNPS, --snps=SNPS  Path to the snps vcf file

###Using adjust_individualized_annotations.py
    ./adjust_individualized_annotations.py --help
    Usage: adjust_individualized_annotations.py [options] annotation-file DO-identifier L-or-R

    Options:
      --version             show program's version number and exit
      -h, --help            show this help message and exit
      -c CHRCOL, --chr-column=CHRCOL
                            Column with chromosome info. (Required, one-based.)
      -C COMMENTS, --comments-file=COMMENTS
                            File to write adjustment comments to. (Default: no
                            comments file.)
      -d DIR, --dir=DIR     path to directory containing offset files (default:
                            current directory)
      -e ENDCOL, --end-column=ENDCOL
                            The column containing the position of the end of the
                            feature. (Required, one-based.)
      -n NEAR, --near=NEAR  Report if a feature position is within near bases of
                            an indel. (Optional, default is no report.)
      -o OFN, --output-file=OFN
                            Output filename (Optional, default: stdout)
      -p POSCOLS, --position-columns=POSCOLS
                            Comma-separated list of other columns with positions
                            to be adjusted. (Optional, one-based.)
      -s STARTCOL, --start-column=STARTCOL
                            The column containing the position of the end of the
                            feature. (Required, one-based.)
      -t TYPECOL, --type-column=TYPECOL
                            Column with feature type info. (Optional, one-based.)
      -x EXTENTS, --extents=EXTENTS
                            File containing the strain extents.

##License

This software is licensed under GPL V3 or later.
