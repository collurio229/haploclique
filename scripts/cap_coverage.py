#!/usr/bin/env python3
"""
Caps a BAM file on a specific coverage threshold and masks all non-GpC positions for epigenetics analysis.

Usage: 
  cap_coverage.py [options] [--] <input> [<output>]
  cap_coverage.py -h | --help



Options:
  -h --help                     Show this message
  -c <num> --coverage <num>     Coverage threshold [default: 300]
"""

import pysam
import random
import sys
from docopt import docopt

def rewrite(read):
    q = read.query_qualities
    p = ''
    i = 0

    while(i < len(read.query_sequence)):
        if i+1 < len(read.query_sequence) and read.query_sequence[i] == 'G':
            if  read.query_sequence[i+1] == 'C':
                p += 'GC'
            else:
                p += 'AA'
            i += 2

        elif i+1 < len(read.query_sequence) and read.query_sequence[i] == 'C':
            p += 'AA'
            i += 2
        else:
            p += 'A'
            i += 1

    read.query_qualities = q
    read.query_sequence = p

def main(argv):

    args = docopt(__doc__, argv=argv)

    if not args['<output>']:
        args['<output>'] = str(args['--coverage']) + '_' + args['<input>']

    cov = int(args['--coverage'])

    infile = pysam.AlignmentFile(args['<input>'], 'rb')
    outfile = pysam.AlignmentFile(args['<output>'], 'wb', template=infile)

    reads = []
    coverage_monitor = [0 for i in range(0, infile.lengths[0])]

    for read in infile.fetch():
        reads.append(read)

    random.shuffle(reads)

    for read in reads:
        
        denied = False

        for i in range(read.reference_start, read.reference_end):
            if coverage_monitor[i] > cov:
                denied = True
                break

        if not denied:
            for i in range(read.reference_start, read.reference_end):
                coverage_monitor[i] += 1

            #rewrite(read)

            outfile.write(read)

if __name__ == '__main__':
    main(sys.argv[1:])
