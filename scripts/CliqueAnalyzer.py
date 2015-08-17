#!/usr/bin/env python3
"""
Usage:
  CliqueAnalyzer.py <logfile>
"""

from docopt import docopt
import re
import sys

def parseLogFile(logfile):
    all_edges = [dict()]
    all_cliques = [dict()]

    ct = 0
    edge = True
    with open(logfile, 'r') as fe:
        for line in fe:
            if (edge):
                m = re.match('---', line)
            else:
                m = re.match(r'>--(\d+)--<', line)

            if m:
                if(edge):
                    edge = False
                    all_cliques.append(dict())
                else:
                    ct = int(m.group(1))
                    all_edges.append(dict())
                    edge = True

                continue

            if edge:
                n = re.match(r'(\d+) -> (.*)', line)
                if n:
                    i = int(n.group(1))
                    all_edges[ct][i] = set([int(s) for s in re.split(' ', n.group(2))])
            else:
                n = re.match(r'(\d+): (.*)', line)
                if n:
                    i = int(n.group(1))
                    all_cliques[ct][i] = set([int(s) for s in re.split(' ', n.group(2))])
    
    return all_edges, all_cliques

def verifyCliques(edges, cliques):
    fail = False
    
    includeall = set()
    for v in cliques.values():
        includeall |= v

    bignode = 0
    for k, v in cliques.items():
        s = v.copy()

        for node in v:
            if node > bignode:
                bignode = node

            if node in edges:
                s &= (edges[node] | {node})
            else:
                s &= {node}
        if s < v:
            print(k, "is not a Clique:", v-s, "aren't included in all nodes")
            fail = True
        elif s > v:
            print(k, "was not maximal:", s-v, "could be added")
            fail = True
#        else:
#            print(k, "was ok")

    allnodes = set(range(0, bignode + 1))
    if includeall < allnodes:
        print(allnodes, includeall, "aren't included in any clique!")
        fail = True
    if fail:
        print("Not all cliques passed!")
        sys.exit(1)
    else:
        print("All cliques passed!")

def main(argv):
    args = docopt(__doc__, argv=argv)

    all_edges, all_cliques = parseLogFile(args['<logfile>'])

    for i in range(0, len(all_edges)):
        verifyCliques(all_edges[i], all_cliques[i])

if __name__ == '__main__':
    main(sys.argv[1:])
