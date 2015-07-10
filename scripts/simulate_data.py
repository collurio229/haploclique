#!/usr/bin/env python3
"""
This script generates simulated haplotype data and Illumina reads from a
reference DNA sequence.

Usage: 
  simulate_data.py (--snp | --ins | --del ) [options] [--] <input> [<output>]
  simulate_data.py -h | --help

  --snp     Generate new haplotypes through substituting random nucleotides
  --ins     Generate new haplotypes through inserting random sequences
  --del     Generate new haplotypes through deletion of a random sequence part

  <input>   The reference sequence in FASTA format
  <output>  The destination of the packed data as tar.gz archive.
            The archive will include a bam file with simulated Illumina reads,
            the reference and all haplotype sequences in FASTA format and
            a log file with all used options.

Options:
  -h --help                     Show this message
  -m <num>, --mean <num>        mean value for length of indels or number of 
                                nucleotide substitutions [default: 1.0]
  -s <num>, --sigma <num>       standard deviation for length of indels or
                                number of nucleotide substitions [default: 0.0]
  -n <num>, --number <num>      number of created haplotypes [default: 1]
  -c <num>, --coverage <num>    average coverage for simulated reads
                                [default: 30]
  -r <n>, --read_length <n>     average read length [default: 100]
  --seed <num>                  seed for random number generator
  --no_ref                      don't include reference sequence in
                                read simulation
  --validate                    validate created sam file with picard
  --error                       print errors of subprocesses to STDERR
  --tmp <dir>                   Write temporary files to dir
"""

from subprocess import call
from docopt import docopt
import random
import tarfile
import tempfile
import sys
import os

try:
    from subprocess import DEVNULL
except ImportError:
    DEVNULL = open(os.devnull, 'wb')

from subprocess import STDOUT
from progressdone import *
from Sequence import *

class ExternalError(Exception):
    """Raise this error if an external process fails.

    This class is a standard implementation of an error class 
    and should be used if a program called with subprocess fails,
    e.g. return value is not 0"""

    def __init__(self, message):
        self.message = message

    def __str__(self):
        return repr(self.message)

def main(argv):
    """This script generates simulated haplotype data and Illumina reads from a reference genome.
    
    Use ['-h'] or ['--help--'] as arguments to display a usage message
    or look into the module __doc__"""

    args = docopt(__doc__, argv=argv)

    if args['--seed']:
        random.seed(args['--seed'])

    if args['--error']:
        STD = None        
    else:
        STD = DEVNULL

    if args['--tmp']:
        TMPDIR = args['--tmp']
    else:
        TMPDIR = '/tmp'

    # convert parsed arguments from string to int or float
    args['--mean'] = float(args['--mean'])
    args['--sigma'] = float(args['--sigma'])
    args['--number'] = int(args['--number'])
    args['--coverage'] = int(args['--coverage'])
    args['--read_length'] = int(args['--read_length'])

    progress('Reading input sequence', False)

    ref_genome = readFASTA(args['<input>'])[0]

    base = os.path.basename(args['<input>'])
    filename = os.path.splitext(base)[0]

    progress('Creating output file', False)

    if not args['<output>']:
        args['<output>'] = 'sht_' + filename + '.tar.gz'

    tar = tarfile.open(args['<output>'], 'w:gz')
    tar.add(args['<input>'], arcname='ref_' + base)

    progress('Writing logfile', False)

    with tempfile.NamedTemporaryFile(prefix= 'log_', suffix='.log', delete=False, dir=TMPDIR) as logfile:
        log = logfile.name        

        if args['--snp']:
            logfile.write(bytes('method: snp\n', 'UTF-8'))
        elif args['--ins']:
            logfile.write(bytes('method: ins\n', 'UTF-8'))
        elif args['--del']:
            logfile.write(bytes('method: del\n', 'UTF-8'))

        logfile.write(bytes('mean: ' + str(args['--mean']) + '\n', 'UTF-8'))
        logfile.write(bytes('sigma: ' + str(args['--sigma']) + '\n', 'UTF-8'))
        logfile.write(bytes('number: ' + str(args['--number']) + '\n', 'UTF-8'))
        logfile.write(bytes('seed: ' + str(args['--seed']) + '\n', 'UTF-8'))
        logfile.write(bytes('coverage: ' + str(args['--coverage']) + '\n', 'UTF-8'))
        
        if args['--no_ref']:
            logfile.write(bytes('no_ref\n', 'UTF-8'))

        if args['--validate']:
            logfile.write(bytes('validate\n', 'UTF-8'))

    tar.add(log, arcname= 'args_' + filename + '.log')

    os.unlink(log)

    progress('Creating new haplotypes')
    subprogress_init_counter()

    n = args['--number']

    haplofiles = []

    for i in range(0, n):

        subprogress(n, 'Creating haplotype %%')

        seq = ref_genome.replicate()

        mean = args['--mean']
        sigma = args['--sigma']

        if mean < 1:
            mean *= len(ref_genome.sequence)
            sigma *= len(ref_genome.sequence)

        if args['--snp']:
            for j in range(0, int(random.gauss(mean, sigma))):

                seq.create_snp(dist_fn=lambda self: random.choice(['A', 'T', 'G', 'C']))

        elif args['--ins']:
            seq.create_insertion(mean, sigma)

        elif args['--del']:
            seq.create_deletion(mean, sigma)

        f = tempfile.NamedTemporaryFile(suffix='.fasta', delete=False, dir=TMPDIR)
        seq.writeFASTA(f)
        haplofiles.append(f.name)
        f.close()
        
        tar.add(haplofiles[-1], arcname='ht' + str(i) + '_' + filename + '.fasta')

    done('Creating new haplotypes')

    if not args['--no_ref']:
        haplofiles.append(args['<input>'])
        n += 1

    samfiles = []

    progress('Simulating sequence reads')
    subprogress_init_counter()

    subprogress(4*n+1, 'Creating bwa index')

    if call(['bwa', 'index', '-p', 'tmp_index', args['<input>']], stdout=STD, stderr=STD):
        raise ExternalError('Failed to create an index for alignment')

    ct = 0
    for hf in haplofiles:

        subprogress(4*n+1, 'Simulating reads for haplotype %%', fn_filter=lambda x: x//4)

        read1 = tempfile.NamedTemporaryFile(prefix='per1_', suffix='.fastq', delete=False, dir=TMPDIR)
        read2 = tempfile.NamedTemporaryFile(prefix='per2_', suffix='.fastq', delete=False, dir=TMPDIR)
        raw_sam = tempfile.NamedTemporaryFile(prefix='sim_', suffix='.sam', delete=False, dir=TMPDIR)
        end_sam = tempfile.NamedTemporaryFile(prefix='aligned_', suffix='.sam', delete=False, dir=TMPDIR)

        samfiles.append(end_sam.name)
        raw_sam.close()
        read1.close()
        read2.close()

        if call(['java', '-jar', '-Xmx2048m', '-Djava.io.tmpdir=' + TMPDIR, 
                '../bin/SimSeq.jar',
                '-1', str(args['--read_length']),
                '-2', str(args['--read_length']),
                '--out', raw_sam.name,
                '--reference', hf,
                '--read_number', str(int(args['--coverage'] * len(ref_genome.sequence) / (2*args['--read_length']) )),
                '--read_prefix', 'SimSeq_ht' + str(ct) + '_',
                '--error', '../data/miseq_250bp.txt'], stdout=STD, stderr=STD):
            raise ExternalError('SimSeq failed for ' + hf)
    
        subprogress(4*n+1, 'Creating sequence dictionary for haplotype %%', fn_filter=lambda x: x//4)

        if call(['java', '-jar', '-Xmx2g', '../bin/picard.jar', 'CreateSequenceDictionary', 'REFERENCE=' + hf, 'OUTPUT=' + TMPDIR + '/dict.sam', 'VERBOSITY=ERROR'], stdout=STD, stderr=STD):
            raise ExternalError('Picard failed to create a sequence dictionary for ' + hf)

        with open(TMPDIR + '/dict.sam', 'a') as header, open(raw_sam.name, 'r') as sam_input:
            for line in sam_input:
                header.write(line)

        subprogress(4*n+1, 'Converting reads of haplotype %% to fastq format', fn_filter=lambda x: x//4)

        if call(['java', '-jar', '-Xmx2g',
                '../bin/picard.jar', 'SamToFastq',
                'INPUT=' + TMPDIR + '/dict.sam',
                'FASTQ=' + read1.name,
                'SECOND_END_FASTQ=' + read2.name,
                'VALIDATION_STRINGENCY=LENIENT',
                'VERBOSITY=ERROR'], stdout=STD, stderr=STD):
            raise ExternalError('Picard failed to convert ' + raw_sam.name + ' to fastq')

        subprogress(4*n+1, 'Aligning reads of haplotype %% to reference genome', fn_filter=lambda x: x//4-1)

        if call(['bwa', 'mem', 'tmp_index', read1.name, read2.name], stdout=end_sam, stderr=STD):
            raise ExternalError('bwa failed to align ' + read1.name + ', ' + read2.name)

        end_sam.close()

        os.unlink(raw_sam.name)
        os.unlink(read1.name)
        os.unlink(read2.name)
        os.unlink(TMPDIR + '/dict.sam')

        if hf != args['<input>']:
            os.unlink(hf)

        ct += 1

    os.unlink('tmp_index.amb')
    os.unlink('tmp_index.ann')
    os.unlink('tmp_index.bwt')
    os.unlink('tmp_index.pac')
    os.unlink('tmp_index.sa')

    done('Simulating sequence reads')

    progress('Merging sam files')
    subprogress_init_counter()

    samf = tempfile.NamedTemporaryFile(prefix='merge_', suffix='.sam', delete=False, dir=TMPDIR)
    samf.close()

    options = ['java', '-jar', '-Xmx6g', '../bin/picard.jar', 'MergeSamFiles', 'OUTPUT=' + samf.name, 'VERBOSITY=ERROR', 'TMP_DIR=' + TMPDIR]

    for sf in samfiles:
        options.append('INPUT=' + sf)

    if call(options, stdout=STD, stderr=STD):
        raise ExternalError('Picard failed at merging all sam files')

    done('Merging sam files')

    for sf in samfiles:
        os.unlink(sf)

    if args['--validate']:
        progress('Validating merged sam file')

        if call(['java', '-jar', '-Xmx2g', '../bin/picard.jar','ValidateSamFile', 'INPUT=' + samf.name]):
            raise ExternalError('Validation of merged sam file failed')

        done('Validating merged sam file')

    progress('Converting reads to bam format', False)

    bamf = tempfile.NamedTemporaryFile(prefix='res_', suffix='.bam', delete=False, dir=TMPDIR)
    bamf.close()

    options = ['java', '-jar', '-Xmx6g', '../bin/picard.jar', 'SamFormatConverter', 'INPUT=' + samf.name, 'OUTPUT=' + bamf.name, 'VERBOSITY=ERROR', 'TMP_DIR=' + TMPDIR]

    if call(options, stdout=STD, stderr=STD):
        raise ExternalError('Picard failed to convert ' + samf.name + ' to bam')

    progress('Adding read alignment file to output file', False)
    tar.add(bamf.name, arcname='reads_' + filename + '.bam')

    os.unlink(samf.name)
    os.unlink(bamf.name)
    tar.close()

    all_done()

if __name__ == '__main__':
    main(sys.argv[1:])
