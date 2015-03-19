#!/usr/bin/env python3

import tarfile
import argparse
import sys
import re

#TODO Generate and read log file for cmd args in tar.gz archives

def main(argv):
    if len(argv != 1):
        print('usage: run_data.py data.tar.gz')

    archive = tarfile.open(argv[0], 'r:gz')

    for member in archive.getmembers():

        if re.match(r'ref_(.*)\.fasta', member.name):
            ref = archive.extractfile(member)
        elif re.match(r'reads_(.*)\.bam', member.name):
            reads = archive.extractfile(member)
        elif re.match(r'ht\d_(.*)\.fasta', member.name):
            haplofiles.append(archive.extractfile(member))
        elif re.match(r'args_(.*)\.log', member.name):
            log = archive.extractfile(member)
        else:
            raise
if __name__ == '__main__':
    main(sys.argv[1:])
