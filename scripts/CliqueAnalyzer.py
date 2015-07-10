#!/usr/bin/env python3
"""
Usage:
  CliqueAnalyzer.py <edgefile> <cliquefile>
  CliqueAnalyzer.py comp cliques <cl1_file> <cl2_file>
  CliqueAnalyzer.py comp edges <edge1_file> <edge2_file>
"""

from docopt import docopt
import re
import sys

def parseEdgeFile(edgefile):
    nodes = dict()
    edges = dict()

    with open(edgefile, 'r') as fe:
        for line in fe:
            m = re.match(r'(?P<id>\d+):(?P<name>.*) -> (?P<edges>.*)', line)
            if m:
                nodes[m.group('name')] = int(m.group('id'))
                edges[int(m.group('id'))] = set([int(s) for s in re.split(' ', m.group('edges'))])

    return nodes, edges

def parseCliqueFile(cliquefile, nodes=None):
    cliques = dict()

    with open(cliquefile, 'r') as fc:
        for line in fc:
            l = re.split(r'[\s,]', line)

            if len(l) >= 2:
                if nodes != None:
                    l[1:-1] = [nodes[i] for i in l[1:-1]]
                cliques[l[0]] = set(l[1:-1])

    return cliques

def verifyCliques(edges, cliques):
    fail = False

    for k, v in cliques.items():
        s = v.copy()
        for node in v:
            s &= (edges[node] | {node})

        if s < v:
            print(k, "is not a Clique:", v-s, "aren't included in all nodes")
            fail = True
        elif s > v:
            print(k, "was not maximal:", s-v, "could be added")
            fail = True
#        else:
#            print(k, "was ok")

    if fail:
        print("Not all cliques passed!")
        sys.exit(1)
    else:
        print("All cliques passed!")

def compareCliques(clique1, clique2):
    fail = False

    for k, v in clique1.items():
        if v != clique2[k]:
            print(k, ": exclusive reads")
            print("0:", v - clique2[k])
            print("1:", clique2[k] -v)
            fail = True

    if fail:
        sys.exit(1)
    else:
        print("cliques inscribed in given files were equal")

def compareGraphs(edge1, node1, edge2, node2):
    fail = False

    for k, v in edge1.items():
        if v != edge2[k]:
            print(v, "is not equal to", edge2[k], "at", k)
            fail = True

    for k, v in node1.items():
        if v != node2[k]:
            print(v, "is not equal to", node2[k], "at", k)
            fail = True
    if fail:
        sys.exit(1)
    else:
        print("edges inscribed in given files were equal")

def main(argv):
    args = docopt(__doc__, argv=argv)

    if args['comp']:
        if args['cliques']:
            clique1 = parseCliqueFile(args['<cl1_file>'])
            clique2 = parseCliqueFile(args['<cl2_file>'])
            compareCliques(clique1, clique2)
        elif args['edges']:
            edge1, node1 = parseEdgeFile(args['<edge1_file>'])
            edge2, node2 = parseEdgeFile(args['<edge2_file>'])
            compareGraphs(edge1, node1, edge2, node2)
        else:
            print("Something went wrong")
            sys.exit(1)
    else:
        nodes, edges = parseEdgeFile(args['<edgefile>'])
        cliques = parseCliqueFile(args['<cliquefile>'], nodes)
        verifyCliques(edges, cliques)
if __name__ == '__main__':
    main(sys.argv[1:])
