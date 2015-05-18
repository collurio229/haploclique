#!/usr/bin/env python3
"""Run haploclique on paired-end reads in bam format.

This is a simplified python version of the original haploclique-assembly
shell script. Because of that, it lacks some of the sophisticated options
of the original script.

Usage:
  haploclique.py [clever] [options] [--] <reference> <input>
  haploclique.py bronkerbosch [options] [--] <reference> <input>
  haploclique.py -h | --help

  <reference>   Reference genome in fasta format
  <input>       Read alignment in bam format
  clever        Use the original CLEVER clique finder
  bronkerbosch  Use the BronKerbosch algorithm for clique finding

Options:
  -i <num>, --iterations <num>  number of haploclique iterations. [default: 1]
  --no-singletons               don't use singletons for new iterations
  --bin <path>                  path to haploclique binaries
  --cleanup                     cleanup directory after last call
"""

from subprocess import check_call, Popen, PIPE, CalledProcessError
from docopt import docopt
import tempfile
import os
import sys
import re
import shutil
try:
    from subprocess import DEVNULL
except ImportError:
    DEVNULL = open(os.devnull, 'wb')

class ParsingError(Exception):
    """Use this error class, if an error while parsing a (FASTA) file occured"""

    def __init__(self, message):
        self.message = message

    def __str__(self):
        return repr(self.message)

def initialize(reference, alignment, path=''):
    """Generate paired alignment priors from bam file"""

    outf = open("alignment.prior", 'wb')
    
    check_call([path + 'bam-to-alignment-priors', '--ignore_xa', reference, alignment], stdout=outf, stderr=DEVNULL)

    outf.close()

    return outf.name

def haploclique(prior, options, path=''):
    """Calls the haploclique program on file prior.

    prior is the filename of the prior file given to haploclique, which can be generated by initialize
    or by assemble, if you want to input the superreads of a previous haploclique run."""

    ret_vals = (0, 0, 0)

    sortprior(prior)

    with open(prior, 'r') as inf, tempfile.SpooledTemporaryFile(mode='w+') as outf:
        
        try:
            check_call([path + 'haploclique'] + options, stdin=inf, stderr=outf)
        except CalledProcessError:
            outf.seek(0)
            for line in outf:
                print(line)
            
            sys.exit(1)

        outf.seek(0)

        for line in outf:
            m = re.search(r'(?P<cliques>\d+)\/(?P<uniques>\d+)\/(?P<cputime>\d+)', line)

            if m:
                ret_vals = (int(m.group('cliques')), int(m.group('uniques')), int(m.group('cputime')))
 
    check_call([path + 'seeksingletons'])

    return ret_vals 

def assemble(reference, read1, read2, read_single, prior, path=''):
    """Assemble the data from a haploclique run.

    This writes the generated superreads into a sorted prior file, which can be given to 
    another haploclique run."""

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

                check_call([path + 'bam-to-alignment-priors', '--ignore_xa', reference, 'paired.bam'], stdout=priorfile, stderr=DEVNULL)

                os.unlink('paired.bam')

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

                check_call([path + 'bam-to-alignment-priors', '--ignore_xa', '--unsorted', '--single-end', reference, 'single.bam'], stdout=priorfile, stderr=DEVNULL)

                os.unlink('single.bam')

    if (os.path.getsize(read_single) == 0 and os.path.getsize(read1) == 0 and os.path.getsize(read2) == 0):
        raise ParsingError('All clique data files are empty!')

def sortprior(prior):
    """Sort prior file as expected from haploclique"""

    with tempfile.NamedTemporaryFile() as tf:

        shutil.copy2(prior, tf.name)
    
        with open(prior, 'w') as priorfile:
            check_call(['sort', '-k6,6', '-g'], stdin=tf, stdout=priorfile)

def main(argv):
    """Run haploclique on paired-end reads in bam format.

    This script is basically a python wrapping for haploclique and is a simple version
    of the haploclique-assembly shell script.

    Use ['-h'] or ['--help'] as arguments to display a usage message
    or look into the module __doc__"""

    args = docopt(__doc__, argv=argv)

    if args['--bin']:
        path = args['--bin']
        if path[-1] != '/':
            path += '/'
    else:
        path = ''

    if args['bronkerbosch']:
        algorithm = 'bronkerbosch'
    else:
        algorithm = 'clever'

    prior = initialize(args['<reference>'], args['<input>'], path)
    
    for i in range(0, int(args['--iterations'])):

        try:
            os.unlink('data_cliques_paired_R1.fastq')
            os.unlink('data_cliques_paired_R2.fastq')
            os.unlink('data_cliques_single.fastq')
            os.unlink('data_clique_to_reads.tsv')
        except OSError:
            pass

        if not args['--no-singletons'] and os.path.exists('singles.prior'):
            check_call(r'sort -u -t" " -k1,1 singles.prior >> ' + prior, shell=True, stderr=DEVNULL)
            os.unlink('singles.prior')

        cliques, uniques, cputime = haploclique(prior, [algorithm, '-L', '100'], path)

        print(cliques, uniques, cputime)

        if cliques == 0 or uniques < 5:
            print('Converged after', i+1, 'runs')
            break

        assemble(args['<reference>'], 'data_cliques_paired_R1.fastq', 'data_cliques_paired_R2.fastq', 'data_cliques_single.fastq', prior, path)

    # Write the found haplotypes into a FASTA file
    check_call(path + r"alignmentparser | sort -r -n -k 2 | sed 's/x /_/g;' | tr ' ' '\n\r' > quasispecies.fasta", shell=True)

    if args['--cleanup']:
        os.unlink('data_cliques_paired_R1.fastq')
        os.unlink('data_cliques_paired_R2.fastq')
        os.unlink('data_cliques_single.fastq')
        os.unlink('data_clique_to_reads.tsv')
        os.unlink(prior)
        if os.path.exists('singles.prior'):
            os.unlink('singles.prior')

if __name__ == '__main__':
    main(sys.argv[1:])
