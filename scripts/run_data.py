#!/usr/bin/env python3

from subprocess import call
import tarfile
import argparse
import sys
import re
import tempfile
import os

from progressdone import *
from Sequence import *

class UnknownFileError(Exception):
    def __init__(self, message):
        self.message = message

    def __str__(self):
        return repr(self.message)

class ExternalError(Exception):
    def __init__(self, message):
        self.message = message

    def __str__(self):
        return repr(self.message)

def main(argv):
    parser = argparse.ArgumentParser(description='Run and compare data sets on haploclique. The tar.gz has to include the following files with specific prefixes:\nReference sequence: ref_name.fasta\nReads in bam format: reads_name.bam\nHaplotype sequences: ht[num]_name.fasta\n')
    parser.add_argument('input', help='archive including reference genome, reads and included haplotypes.', metavar='data_archive.tar.gz')
    parser.add_argument('--path_var', help='Set and unset the path variable to evade installing haploclique', default=False, action='store_true')
    parser.add_argument('--iterations', help='number of haploclique iterations. Default: 3', default=3)

    args = parser.parse_args(argv)

    archive = tarfile.open(args.input, 'r:gz')

    path = tempfile.mkdtemp(prefix='hcl_')

    haplofiles = []

    for member in archive.getmembers():

        if re.match(r'ref_(.*)\.fasta', member.name):
            archive.extract(member, path)
            ref = member.name
        elif re.match(r'reads_(.*)\.bam', member.name):
            archive.extract(member, path)
            reads = member.name
        elif re.match(r'ht(\d+)_(.*)\.fasta', member.name):
            archive.extract(member, path)
            haplofiles.append(member.name)
        elif re.match(r'args_(.*)\.log', member.name):
            log = archive.extractfile(member)
        else:
            raise UnknownFileError('Unrecognized file. Check prefix or file extensions of included files')

    # set the PATH variable, so haploclique has not to be installed
    if args.path_var:
        pathvar = os.environ['PATH']
        os.environ['PATH'] = '/home/stud/lanber/usr/lib/jvm/java-7-openjdk-amd64/jre/bin:/home/stud/lanber/usr/bin:' + pathvar + ':' + os.path.realpath('../bin') + ':' + os.path.realpath('.')
        os.environ['SAF'] = os.path.realpath('../bin')

    if call(['samtools', 'index', path + '/' + reads]):
        raise ExternalError('Index creation failed')

    if call(['./haploclique-assembly', '-r', path + '/' + ref, '-i', path + '/' + reads, '-t', str(args.iterations)]):
        raise ExternalError('Haploclique failed')

    quasispecies = readFASTA('./quasispecies.fasta')
    ref_sequences = []

    for hf in haplofiles:
        ref_sequences.extend(readFASTA(path + '/' + hf))

    print('forward error:', compute_best_match_score(ref_sequences, quasispecies))
    print('backward error:', compute_best_match_score(quasispecies, ref_sequences))

    # clean PATH variable
    if args.path_var:
        os.environ['PATH'] = pathvar

if __name__ == '__main__':
    main(sys.argv[1:])
