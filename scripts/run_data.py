#!/usr/bin/env python3
"""This script runs haploclique on the given reference sequence and
bam alignment and compares the results with the given haplotype sequences.

Usage:
  run_data.py [options]

Options:
  -p <path>, --path <path>      path to haploclique binaries [default: ../bin]
"""

from subprocess import call, check_call, CalledProcessError
from docopt import docopt
from shutil import rmtree as remove_directory
import tarfile
import argparse
import sys
import re
import tempfile
import os
import resource
import signal

from progressdone import *
from Sequence import *

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

def write_data(bk, cl, coverage, filename):
    with open(filename, 'w') as data:
        c = ''
        for i in coverage:
            c += ' ' + str(i)

        print(r'bronkerbosch:' + c, file=data)
        
        for line in bk:
            print(line, file=data)

        print(r'clever:' + c, file=data)
        
        for line in cl:
            print(line, file=data)

def main(argv):
    """This script runs haploclique on the given reference sequence and
    bam alignment and compares the results with the given haplotype sequences.

    Use ['-h'] or ['--help--'] as arguments to display a usage message
    or look into the module __doc__"""

    args = docopt(__doc__, argv=argv)

    ARCHIVES = ['sht_arabis_short', 'sht_arabis_large', 'sht_HIV']
    COVERAGES = [450]#[75, 150, 225, 300, 450]
    snp = '_01'

    bk = []
    cl = []

    for archive in ARCHIVES:

        bk.append(archive)
        cl.append(archive)

    ct = 0

    for c in COVERAGES:
        for archive in ARCHIVES:
            print('Processing', archive + '_' + str(c) + '.tar.gz')

            try:
                (cl_time, cl_cliques) = execute(archive + '_' + str(c) + snp + '.tar.gz', args, False)
            except KeyboardInterrupt:
                write_data(bk, cl, COVERAGES, 'data.tex')
                sys.exit(1)
            except CalledProcessError as e:
                #write_data(bk, cl, COVERAGES, 'data.tex')
                (cl_time, cl_cliques) = ('nan', 0)
                print(e)
                #sys.exit(1)

            try:
                (bk_time, bk_cliques) = execute(archive + '_' + str(c) + snp + '.tar.gz', args, True)
            except KeyboardInterrupt:
                write_data(bk, cl, COVERAGES, 'data.tex')
                sys.exit(1)
            except CalledProcessError as e:
                #write_data(bk, cl, COVERAGES, 'data.tex')
                (bk_time, bk_cliques) = ('nan', 0)
                print(e)
                #sys.exit(1)

#                if (bk_cliques != cl_cliques):
#                    print('bk:')
#                    for count in bk_cliques:
#                        print(count)
#                    print('cl:')
#                    for count in cl_cliques:
#                        print(count)
#                    sys.exit(1)

            bk[ct] += ' ' + str(bk_time)
            cl[ct] += ' ' + str(cl_time)

            ct += 1

        ct = 0

    write_data(bk, cl, COVERAGES, 'data.tex')

def execute(dataset, args, bk):
    progress('Opening archive', False)

    archive = tarfile.open(dataset, 'r:gz')

    path = tempfile.mkdtemp(prefix='hcl_')

    for member in archive.getmembers():

        if re.match(r'ref_(.*)\.fasta', member.name):
            #archive.extract(member, path)
            ref = member.name
        elif re.match(r'reads_(.*)\.bam', member.name):
            archive.extract(member, path)
            reads = member.name
        elif re.match(r'ht(\d+)_(.*)\.fasta', member.name):
            pass
            #archive.extract(member, path)
            #haplofiles.append(member.name)
        elif re.match(r'args_(.*)\.log', member.name):
            logfile = archive.extractfile(member)
        else:
            raise UnknownFileError('Unrecognized file. Check prefix or file extensions of included files')

    log = parseLog(logfile)

    #if not log['no_ref']:
    #    haplofiles.append(ref)

    # set the PATH variable to local haploclique binary
    pathvar = os.path.realpath(args['--path'])

    progress('Creating bam index', False)

    if call(['samtools', 'index', path + '/' + reads]):
        raise ExternalError('Index creation failed')

    progress('Calling haploclique')

    options = []

    if bk:
        options += ['bronkerbosch']

#    options += ['-s', '0.0', path + '/' + reads]
    options += ['-n', path + '/' + reads]

    time = 0.0

    with tempfile.TemporaryFile() as tmpout:
        check_call([pathvar + '/haploclique'] + options, stdout=tmpout)

        tmpout.seek(0)

        cliques = []

        for line in tmpout:
            m = re.match(r'time: (?P<time>[0-9.]+)', line.decode())
            n = re.match(r'\d+: (?P<clique>\d+)', line.decode())
            if m:
                time = float(m.group('time'))
            elif n:
                cliques.append(int(n.group('clique')))

    print('needed', time, 'seconds')

    done('Calling haploclique')

    return (time, cliques)

if __name__ == '__main__':
    main(sys.argv[1:])
