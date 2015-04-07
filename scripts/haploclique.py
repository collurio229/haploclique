#!/usr/bin/env python3
"""
Run haploclique on paired-end reads in bam format.

Usage:
  haploclique.py [options] [--] <reference> <input>
  haploclique.py -h | --help

  <reference>   Reference genome in fasta format
  <input>       Read alignment in bam format

Options:
  -i <num>, --iterations <num>  number of haploclique iterations. [default: 1]
"""

from subprocess import check_call, DEVNULL
from docopt import docopt
import tempfile
import os

def initialize(reference, alignment):
    """Generate paired alignment priors from bam file"""

    outf = tempfile.NamedTemporaryFile(suffix='.prior', delete=False)
    check_call(['bam-to-alignment-priors', '--ignore_xa', reference, alignment], stdout=outf, stderr=DEVNULL)
    outf.close()

    return outf.name

def haploclique(prior):
    """Calls the haploclique program"""

    with open(prior, 'r') as inf:
        
        check_call(['haploclique', '-L', '100'], stdin=inf)

def assemble(reference, read1, read2, prior):
    check_call(['bwa', 'index', reference], stdout=DEVNULL, stderr=DEVNULL)

    samfile = tempfile.NamedTemporaryFile(suffix='.sam', delete=False)
    check_call(['bwa', 'mem', '-t', '8', reference, read1, read2], stdout=samfile)
    samfile.close()

    check_call(['samtools', 'faidx', reference], stdin=DEVNULL, stderr=DEVNULL)

    bamfile = tempfile.NamedTemporaryFile(suffix='.bam', delete=False)
    check_call(['samtools', 'view', '-q', '1', '-F', '4', '-bt', reference + '.fai', samfile.name], stdout=bamfile, stderr=DEVNULL)
    bamfile.close()

    check_call(['samtools', 'sort', bamfile.name, 'reads'], stdin=DEVNULL, stderr=DEVNULL)

    os.unlink(samfile.name)
    os.unlink(bamfile.name)
    with open(prior, 'w') as priorfile:
        check_call(['bam-to-alignment-priors', '--ignore_xa', reference, 'reads.bam'], stdout=priorfile, stderr=DEVNULL)
    
def main(argv):
    """See module doc for usage"""

    args = docopt(__doc__, argv=argv)

    prior = initialize(args['<reference>'], args['<input>'])
    
    haploclique(prior)

    os.unlink(prior)

if __name__ == '__main__':
    main(sys.argv[1:])
