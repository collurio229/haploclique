#!/usr/bin/env python3

from subprocess import call
import random
import re
import argparse
import tarfile
import tempfile
import sys
import os

from progressdone import *

class Sequence:

    def __init__(self, identifier, description, bases=['A', 'G', 'C', 'T']):
        self.identifier = identifier
        self.description = description
        self.sequence = []
        self.bases = bases

        self.base_freqs = dict()
        for b in self.bases:
            self.base_freqs[b] = 0.25

    def replicate(self):
        new_seq = Sequence(self.identifier, self.description)
        new_seq.add_sequence(''.join(self.sequence))

        return new_seq

    def add_sequence(self, seq):
        old_len = len(self.sequence)
        self.sequence.extend(list(seq))

        length = len(self.sequence)

        new_len = len(seq)

        new_freq = dict()
        for key in self.base_freqs.keys():
            new_freq[key] = 0

        for c in seq:
            new_freq[c] += 1
            
        for (base, freq) in new_freq.items():
            self.base_freqs[base] = (old_len * self.base_freqs[base] + freq) / length

    def sample_base(self):
        r = random.random()
        ct = 0

        for (base, freq) in self.base_freqs.items():
            ct += freq
            if r < ct:
                return base

    def create_snp(self, dist_fn=lambda self: self.sample_base(), pos=-1, base='N'):
        if pos < 0 or pos > len(self.sequence):
            pos = random.randint(0, len(self.sequence)-1)

        if base == 'N':
            base = dist_fn(self)

        self.sequence[pos] = base

    def create_insertion(self, mean, std_deviation, pos=-1):
        length = int(random.gauss(mean, std_deviation))

        if pos < 1 or pos > len(self.sequence):
            pos = random.randint(0, len(self.sequence)-1)

        insert = []

        for i in range(0, length):

            insert.append(self.sample_base())

        self.sequence[pos:pos] = insert

    def create_deletion(self, mean, std_deviation, pos=-1):
        length = int(random.gauss(mean, std_deviation))

        if pos < 1 or pos > len(self.sequence)-length:
            pos = random.randint(0, len(self.sequence)-length-1)

        self.sequence[pos:pos+length] = []

    def writeFASTA(self, fileobj):
        fileobj.write(bytes('>' + self.identifier + '|' + self.description + '\n', 'UTF-8'))
        ct = 0
        for base in self.sequence:
            fileobj.write(bytes(base, 'UTF-8'))
            ct += 1
            if ct == 70:
                ct = 0
                fileobj.write(b'\n')

class ParsingError(Exception):
    def __init__(self, message):
        self.message = message

    def __str__(self):
        return repr(self.message)

class ExternalError(Exception):
    def __init__(self, message):
        self.message = message

    def __str__(self):
        return repr(self.message)


def readFASTA(filename):
    with open(filename, 'r') as f:
        header = f.readline()

        m = re.search('>(?P<identifier>\w*)\|(?P<description>.*)', header)
        if not m:
            raise ParsingError('FASTA header was not valid')

        seq = Sequence(m.group('identifier'), m.group('description'))

        for line in f:
            seq.add_sequence(line.rstrip())

    return seq



def main(argv):
    parser = argparse.ArgumentParser(description='Generate haplotypes from a reference genome and create simulated illumina reads with SimSeq')
    parser.add_argument('-i', '--input', help='the reference genome in FASTA format (REQUIRED)', required=True, metavar='reference.fasta')
    parser.add_argument('-o', '--output', help='destination of the generated data tar.gz archive. It contains the reference sequence, the generated haplotypes and the simulated sam file.', metavar='output.tar.gz')
    parser.add_argument('--snp', help='generate new haplotypes through substituting random nucleotides.', default=False, action='store_true')
    parser.add_argument('--insert', help='generate new haplotypes through inserting random sequences.', default=False, action='store_true')
    parser.add_argument('--delete', help='generate new haplotypes through deleting random sequences.', default=False, action='store_true')
    parser.add_argument('-m', '--mean', type=int, default=1, help='mean value for length of indels or number of nucleotide substitutions. Default: 1')
    parser.add_argument('-s', '--sigma', type=int, default=0, help='standard deviation for length of indels or number of nucleotide substitions. Default: 0')
    parser.add_argument('-n', '--number', type=int, default=1, help='number of created haplotypes. Default: 1')
    parser.add_argument('--seed', type=int, help='seed for the Random Number Generator')
    parser.add_argument('-c', '--coverage', type=float, default=30.0, help='Average coverage for simulated reads. Default: 30')
    parser.add_argument('--no_ref', help='Don\'t add reference genome to data simulation.', default=True, action='store_false')
    args = parser.parse_args(argv)

    if args.seed:
        random.seed(args.seed)

    progress('Reading input sequence')

    ref_genome = readFASTA(args.input)

    base = os.path.basename(args.input)
    filename = os.path.splitext(base)[0]

    done('Reading input sequence')

    progress('Creating output file')

    if not args.output:
        args.output = filename + '_sht.tar.gz'

    tar = tarfile.open(args.output, 'w:gz')
    tar.add(args.input, arcname='ref_' + base)

    done('Creating output file')

    progress('Creating new haplotypes')
    subprogress_init_counter()

    n = args.number

    haplofiles = []

    if not args.no_ref:
        haplofiles.append(args.input)

    for i in range(0, n):

        subprogress(n, 'Creating haplotype %%')

        seq = ref_genome.replicate()

        mean = args.mean
        sigma = args.sigma

        if mean < 1:
            mean *= len(ref_genome.sequence)
            sigma *= len(ref_genome.sigma)

        if args.snp:
            for j in range(0, int(random.gauss(mean, sigma))):
                seq.create_snp()

        if args.insert:
            seq.create_insertion(mean, sigma)

        if args.delete:
            seq.create_deletion(mean, sigma)

        f = tempfile.NamedTemporaryFile(suffix='.fasta', delete=False)
        seq.writeFASTA(f)
        haplofiles.append(f.name)
        f.close()
        
        tar.add(haplofiles[-1], arcname='ht' + str(i) + '_' + filename + '.fasta')

    done('Creating new haplotypes')

    samfiles = []

    progress('Simulating sequence reads')
    subprogress_init_counter()

    #subprogress(n+1, 'Creating sequence dictionary')

    #call(['java', '-jar', '../bin/picard.jar', 'CreateSequenceDictionary', 'REFERENCE=' + args.input, 'OUTPUT=dict.sam'])

    for hf in haplofiles:

        subprogress(n+1, 'Simulating reads for haplotype %%')

        read1 = tempfile.NamedTemporaryFile(prefix='per1_', suffix='.fastq', delete=False)
        read2 = tempfile.NamedTemporaryFile(prefix='per2_', suffix='.fastq', delete=False)
        raw_sam = tempfile.NamedTemporaryFile(prefix='sim_', suffix='.sam', delete=False)
        end_sam = tempfile.NamedTemporaryFile(prefix='aligned_', suffix='.sam', delete=False)

        samfiles.append(end_sam.name)
        raw_sam.close()
        end_sam.close()
        read1.close()
        read2.close()

        if call(['java', '-jar', '-Xmx2048m', '../bin/SimSeq.jar',
                '--out', raw_sam.name,
                '--reference', hf,
                '--read_number', str(int(args.coverage * len(ref_genome.sequence))),
                '--error', '../data/miseq_250bp.txt']):
            raise ExternalError('SimSeq failed for ' + hf)
    
        if call(['java', '-jar', '../bin/picard.jar', 'SamToFastq', 'INPUT=' + raw_sam.name, 'FASTQ=' + read1.name, 'SECOND_END_FASTQ=' + read2.name]):
            raise ExternalError('Picard failed to convert ' + raw_sam.name + ' to fastq')

        os.unlink(raw_sam.name)
        os.unlink(read1.name)
        os.unlink(read2.name)

    for hf in haplofiles:
        if hf != args.input:
            os.unlink(hf)

    done('Simulating sequence reads')

    progress('Merging sam files')
    subprogress_init_counter()

    samf = tempfile.NamedTemporaryFile(prefix='merge_', suffix='.sam', delete=False)
    samf.close()

    options = ['java', '-jar', '../bin/picard.jar', 'MergeSamFiles', 'OUTPUT=' + samf.name, 'VALIDATION_STRINGENCY=LENIENT']
    options_sort = ['java', '-jar', '../bin/picard.jar', 'SortSam', 'SORT_ORDER=coordinate', 'VALIDATION_STRINGENCY=LENIENT']

    for sf in samfiles:
        #subprogress(n+1, 'Sorting sam file %%')
        #if call(options_sort + ['INPUT=' + sf, 'OUTPUT=' + sf]):
        #    raise ExternalError('Picard failed at sorting ' + sf)
        options.append('INPUT=' + sf)

    #subprogress(n+1, 'Merging all sam files')

    if call(options):
        raise ExternalError('Picard failed at merging all sam files')

    done('Merging sam files')

    for sf in samfiles:
        os.unlink(sf)

    progress('Converting reads to bam format')

    bamf = tempfile.NamedTemporaryFile(prefix='res_', suffix='.bam', delete=False)
    bamf.close()

    options = ['java', '-jar', '../bin/picard.jar', 'SamFormatConverter', 'INPUT=' + samf.name, 'OUTPUT=' + bamf.name]

    if call(options):
        raise ExternalError('Picard failed to convert ' + samf.name + ' to bam')

    progress('Adding read alignment file to output file')
    tar.add(bamf.name, arcname='reads_' + filename + '.bam')

    done('Adding reads to output file')

    os.unlink(samf.name)
    os.unlink(bamf.name)
    tar.close()

    all_done()

if __name__ == '__main__':
    main(sys.argv[1:])