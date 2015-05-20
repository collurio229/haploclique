#!/usr/bin/env python3
"""This script runs haploclique on the given reference sequence and
bam alignment and compares the results with the given haplotype sequences.

Usage:
  run_data.py [options] <input>

  <input>   tar.gz archive which was either generated with simulate_data or
            follows this naming conventions:
              ref_name.fasta   - reference sequence
              reads_name.bam   - reads as bam alignment
              ht<i>_name.fasta - haplotype sequence number <i>
              args_name.log    - logfile containing used arguments (optional)

Options:
  -i <num>, --iterations <num>  number of haploclique iterations [default: 5]
  -p <path>, --path <path>      path to haploclique binaries [default: ../bin]
  -m --metric                   Compute match metric between quasispecies
                                and haplotypes
  -b --bronkerbosch             Use Bron-Kerbosch algorithm
"""

from subprocess import call
from docopt import docopt
from shutil import rmtree as remove_directory
import tarfile
import argparse
import sys
import re
import tempfile
import os

from progressdone import *
from Sequence import *
import haploclique

class UnknownFileError(Exception):
    """Raise this error if a file doesn't apply to naming conventions.

    This class is a standard implementation of an error class 
    and should be used if a file in the given tar.gz archive doesn't
    follow the naming conventions. For naming conventions use -h"""

    def __init__(self, message):
        self.message = message

    def __str__(self):
        return repr(self.message)

class ExternalError(Exception):
    """Raise this error if an external process fails.

    This class is a standard implementation of an error class 
    and should be used if a program called with subprocess fails,
    e.g. return value is not 0"""

    def __init__(self, message):
        self.message = message

    def __str__(self):
        return repr(self.message)

def parseLog(logfile):
    """Parses a log file into a dictionary.

    The log file can contain information about the method used to
    create artificial haplotypes, the numerical arguments used and
    set binary flags.
    
    This function returns a dictionary containing the given information.
    """
    
    regexes = {'method': re.compile(r'method: (?P<method>\w+)'),
               'numeric': re.compile(r'(?P<name>\w+): (?P<digit>\d+)'),
               'binary': re.compile(r'(?P<name>no_ref|validate)')}

    log = {'snp': False, 'ins': False, 'del': False, 'no_ref': False, 'validate': False}

    for line in logfile:

        for key, regex in regexes.items():
            m = regex.match(str(line))

            if (m):
                if key == 'method':
                    log[m.group('method')] = True
                    break
                elif key == 'numeric':
                    log[m.group('name')] = m.group('digit')
                    break
                else:
                    log[m.group('name')] = True
                    break
                    
                    
    return log

def main(argv):
    """This script runs haploclique on the given reference sequence and
    bam alignment and compares the results with the given haplotype sequences.

    Use ['-h'] or ['--help--'] as arguments to display a usage message
    or look into the module __doc__"""

    args = docopt(__doc__, argv=argv)

    progress('Opening archive', False)

    archive = tarfile.open(args['<input>'], 'r:gz')

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
            logfile = archive.extractfile(member)
        else:
            raise UnknownFileError('Unrecognized file. Check prefix or file extensions of included files')

    log = parseLog(logfile)

    if not log['no_ref']:
        haplofiles.append(ref)

    # set the PATH variable to local haploclique binary
    pathvar = args['--path']
    os.environ['SAF'] = os.path.realpath('../bin')

    progress('Creating bam index', False)

    if call(['samtools', 'index', path + '/' + reads]):
        raise ExternalError('Index creation failed')

    progress('Calling haploclique')

    algorithm = 'clever'

    if args['--bronkerbosch']:
        algorithm = 'bronkerbosch'

    haploclique.main([algorithm, '--iterations', str(args['--iterations']), '--cleanup', '--no-singletons', '--bin', pathvar, '--', path + '/' + ref, path + '/' + reads])

    done('Calling haploclique')

    quasispecies = readFASTA('./quasispecies.fasta')
    ref_sequences = []

    if(args['--metric']):
        progress('Computing match metric', False)

        for hf in haplofiles:
            ref_sequences.extend(readFASTA(path + '/' + hf))

        print('forward error:', compute_best_match_score(ref_sequences, quasispecies))
        print('backward error:', compute_best_match_score(quasispecies, ref_sequences))

    remove_directory(path)

    all_done()

if __name__ == '__main__':
    main(sys.argv[1:])
