#Seqnature

The programs in this repository reflect as-yet unpublished work.
Please contact Dr. Steven Munger, Ph.D. (steven.munger at jax.org) for further
information.

The programs' author is Al Simons (al.simons at jax.org).  Please contact him
to report any problems with the programs.

_Seqnature_ is a set of programs to create custom sequences from a reference 
sequence and a pair of VCF files containing SNP and Indel information.  
It was created for, and most heavily tested on, mouse.

If you have gtf files with annotations based on the same reference, the
programs can create an updated gtf, reflecting the Indels.

##Creating new sequence and annotations

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

##Using build_new_sequence_from_vcfs.py

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
      -r REFERENCE, --reference=REFERENCE
                            Path to the reference file used as the base.
      -s SNPS, --snps=SNPS  Path to the snps vcf file

##Using adjust_annotations.py

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

##License

This software is licensed under GPL V3 or later.