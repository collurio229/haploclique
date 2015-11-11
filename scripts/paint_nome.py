#!/usr/bin/env python3
"""
Caps a BAM file on a specific coverage threshold and masks all non-GpC positions for epigenetics analysis.

Usage: 
  paint_nome.py [options] [--] <ref> <seq>
  paint_nome.py -h | --help



Options:
  -h --help                     Show this message
  -c <num> --coverage <num>     Coverage threshold [default: 300]
"""

import random
import sys
import re
import math
from docopt import docopt

import pygame
from pygame.locals import *
from pygame.gfxdraw import *
import Sequence

def score(a, b):

    convert = {'A': 0, 'T': 1, 'C': 2, 'G': 3, '-': 4, 'N': 4}

#             A  T  C  G  -
    table = [[1, 0, 0, 0, 0], #A
             [0, 1, 0, 0, 0], #T
             [0, 1, 2, 0, 0], #C
             [0, 0, 0, 2, 0], #G
             [0, 0, 0, 0, 0]] #-

    x = convert[a]
    y = convert[b]
    return table[x][y]

def align(ref, seq):
    s = ''
    for c in seq:
        if c != 'N':
            s += c

    seq = s

    l = [[0 for j in range(0, len(ref) + 1)] for i in range(0, len(seq) + 1)]

    for i in range(0, len(seq)):
        for j in range(0, len(ref)):
            x = l[i][j] + score(ref[j], seq[i])
#            y = l[i][j+1] + 0
            z = l[i+1][j]

            l[i+1][j+1] = max(x, z)

    maxx = l[len(seq)][0]
    pos = 0

    for i in range(1, len(ref) + 1):
        if l[len(seq)][i] > maxx:
            maxx = l[len(seq)][i]
            pos = i 

    r = ''
    s = ''

    i = len(seq)
    while(i != 0 and pos != 0):
        x = l[i-1][pos-1] + score(ref[pos-1], seq[i-1])
#        y = l[i-1][pos] + 0
        z = l[i][pos-1]

        if i > 0 and pos > 0 and x == l[i][pos]:
            pos -= 1
            i -= 1
            s += seq[i]            
#        elif i > 0 and y == l[i][pos]:
#            i -= 1
#            r += '-'
#            s += seq[i]
        elif pos > 0 and z == l[i][pos]:
            pos -= 1
            s += '-'
        else:
            print('something went wrong')
            sys.exit(-1)

#    if pos == 0:
#        while( i != 0):
#            i -= 1
#            r += '-'
#            s += seq[i]
#    else:
    while(pos != 0):
        pos -= 1
        r += ref[pos]
        s += '-'

    return s[::-1]

def paint(ref, refname, seqs):
    pygame.init()

    red = (222, 135, 135)
    green = (170, 222, 135)
    blue = (135, 205, 222)
    yellow = (255, 204, 0)
    black = (0, 0, 0)
    white = (255, 255, 255)

    size = width, height = 1200, 50 * (len(seqs) + 1)

    plot = pygame.Surface(size)
    gc = pygame.Surface(size)

    plot.fill(white)
    gc.set_colorkey(black)
    font = pygame.font.SysFont('Arial', 20)

    y = 75
    
    pygame.draw.line(plot, black, (150 , 25), (150, 30), 2)
    for i in range(100, len(ref)+100, 100):
        pygame.draw.line(plot, black, (150 + 2*(i - 100), 25), (150 + 2 * i, 25), 2)    
        pygame.draw.line(plot, black, (150 + 2*i, 25), (150 +2*i, 30), 2)

    pygame.draw.line(plot, blue, (150, y), (150 + 2*len(ref), y), 10)
    for i in range(0, len(ref) -1):
        if ref[i] == 'G' and ref[i+1] == 'C':
            pygame.draw.rect(gc, yellow, (2*i-2 + 150, 65, 4, 20))


    a = 0

    for (p, s) in seqs:
        y += 50
        pos = 0
        gap = True

        #c, a = circle(a, p)
        t = font.render(("%.1f" % (p*100)) + '%', True, (10, 10, 10))
        #plot.blit(c, (25, y - 50))
        plot.blit(t, (50, y-10))
        for i in range(0, len(s)):
            if gap and s[i] == '-':
                continue
            elif not gap and (s[i] == '-' or i == len(s) - 1):
                gap = True
                pygame.draw.rect(plot, green, (2*pos + 150, y - 5, 2*i - 2*pos , 10))
            elif gap and s[i] != '-':
                gap = False
                pos = i

            if i != len(s)-1 and s[i] == 'G' and s[i+1] =='C':
                gap = False
                pygame.draw.rect(gc, yellow, (2*i-2 + 150, y-10, 4, 20))

    screen = pygame.display.set_mode(size)
    screen.blit(plot, (0, 0))
    screen.blit(gc, (0, 0))

    pygame.image.save(screen, refname +'.png') 

    pygame.display.flip()

    while 1:
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                sys.exit()

def main(argv):

    args = docopt(__doc__, argv=argv)

    ref = Sequence.readFASTA(args['<ref>'])[0]
    seqs = Sequence.readFASTA(args['<seq>'])
    
    aligned = []

    for s in seqs:
        seq = align(ref.sequence, s.sequence)
        m = re.search('freq:(0\.\d+)', s.identifier)
        print(seq)
        aligned.append( (float(m.group(1)),seq) )

    paint(ref.sequence, ref.identifier, aligned)

if __name__ == '__main__':
    main(sys.argv[1:])
