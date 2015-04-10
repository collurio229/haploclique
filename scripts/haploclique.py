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

from subprocess import check_call, Popen, PIPE
from docopt import docopt
import tempfile
import os
import sys
import re
try:
    from subprocess import DEVNULL
except ImportError:
    DEVNULL = open(os.devnull, 'wb')

def initialize(reference, alignment):
    """Generate paired alignment priors from bam file"""

    outf = tempfile.NamedTemporaryFile(suffix='.prior', delete=False)
    
    ps = Popen(['../bin/bam-to-alignment-priors', '--ignore_xa', reference, alignment], stdout=PIPE, stderr=DEVNULL)

    check_call(['sort', '-k6,6', '-g'], stdin=ps.stdout, stdout=outf)
    ps.wait()

    outf.close()

    return outf.name

def haploclique(prior):
    """Calls the haploclique program on file prior.

    prior is the filename of the prior file generated by initialize."""

    ret_vals = (0, 0, 0)

    with open(prior, 'r') as inf, tempfile.SpooledTemporaryFile(mode='w+') as outf:
        
        check_call(['../bin/haploclique', '-L', '100'], stdin=inf, stderr=outf)

        outf.seek(0)

        for line in outf:
            m = re.search(r'(?P<cliques>\d+)\/(?P<uniques>\d+)\/(?P<cputime>\d+)', line)

            if m:
                ret_vals = (int(m.group('cliques')), int(m.group('uniques')), int(m.group('cputime')))
 
    return ret_vals 

def assemble(reference, read1, read2, read_single, prior):
    """Assemble the data from a haploclique run."""

    check_call(['bwa', 'index', reference], stdout=DEVNULL, stderr=DEVNULL)
    check_call(['samtools', 'faidx', reference], stdin=DEVNULL, stderr=DEVNULL)
    
    samfile, bamfile = None, None

    with open(prior, 'w') as priorfile:

        if (os.path.exists(read1) and os.path.exists(read2)):
            if (os.path.getsize(read1) != 0 and os.path.getsize(read2) != 0):

                samfile = tempfile.NamedTemporaryFile(suffix='.sam', delete=False)
                try:
                    check_call(['bwa', 'mem', '-t', '8', reference, read1, read2], stdout=samfile, stderr=DEVNULL)
                finally:
                    samfile.close()

                bamfile = tempfile.NamedTemporaryFile(suffix='.bam', delete=False)
                try:
                    check_call(['samtools', 'view', '-q', '1', '-F', '4', '-bt', reference + '.fai', samfile.name], stdout=bamfile, stderr=DEVNULL)
                finally:
                    bamfile.close()

                check_call(['samtools', 'sort', bamfile.name, 'paired'], stdout=DEVNULL, stderr=DEVNULL)

                os.unlink(samfile.name)
                os.unlink(bamfile.name)

                ps = Popen(['../bin/bam-to-alignment-priors', '--ignore_xa', reference, 'paired.bam'], stdout=PIPE, stderr=DEVNULL)

                check_call(['sort', '-k6,6', '-g'], stdin=ps.stdout, stdout=priorfile)
                ps.wait()

        if (os.path.exists(read_single)):
            if (os.path.getsize(read_single) != 0):

                samfile = tempfile.NamedTemporaryFile(suffix='.sam', delete=False)
                check_call(['bwa', 'mem', '-t', '8', reference, read_single], stdout=samfile, stderr=DEVNULL)
                samfile.close()

                bamfile = tempfile.NamedTemporaryFile(suffix='.bam', delete=False)
                check_call(['samtools', 'view', '-q', '1', '-F', '4', '-bt', reference + '.fai', samfile.name], stdout=bamfile)
                bamfile.close()

                check_call(['samtools', 'sort', bamfile.name, 'single'], stdout=DEVNULL, stderr=DEVNULL)

                os.unlink(samfile.name)
                os.unlink(bamfile.name)

                ps = Popen(['../bin/bam-to-alignment-priors', '--ignore_xa', '--unsorted', '--single-end', reference, 'single.bam'], stdout=PIPE, stderr=DEVNULL)

                check_call(['sort', '-k6,6', '-g'], stdin=ps.stdout, stdout=priorfile)
                ps.wait()

                os.unlink('single.bam')

def main(argv):
    """See module doc for usage"""

    args = docopt(__doc__, argv=argv)

    prior = initialize(args['<reference>'], args['<input>'])
    
    for i in range(0, int(args['--iterations'])):

        try:
            os.unlink('data_cliques_paired_R1.fastq')
            os.unlink('data_cliques_paired_R2.fastq')
            os.unlink('data_cliques_single.fastq')
            os.unlink('data_clique_to_reads.tsv')
        except OSError:
            pass

        cliques, uniques, cputime = haploclique(prior)

        if cliques == 0:
            print('Converged after', i+1, 'runs')
            break

        assemble(args['<reference>'], 'data_cliques_paired_R1.fastq', 'data_cliques_paired_R2.fastq', 'data_cliques_single.fastq', prior)

    check_call(r"../bin/alignmentparser | sort -r -n -k 2 | sed 's/x /_/g;' | tr ' ' '\n\r' > quasispecies.fasta", shell=True)

    os.unlink('data_cliques_paired_R1.fastq')
    os.unlink('data_cliques_paired_R2.fastq')
    os.unlink('data_cliques_single.fastq')
    os.unlink('data_clique_to_reads.tsv')
    os.unlink(prior)

if __name__ == '__main__':
    main(sys.argv[1:])
