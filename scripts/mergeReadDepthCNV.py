#! /usr/bin/env python

import sys

cnvkit = sys.argv[1]
cnvpytor = sys.argv[2]
freec = sys.argv[3]
mops = sys.argv[4]
outputFile = sys.argv[5]

def readFile(infile):
    cnv = []
    with open(infile, 'r') as f:
        for x in f:
            if x.startswith('chromosome'):
                continue
            x = x.strip().split('\t')
            if len(x) >= 5:
                x.pop(4)
                cnv.append(x)
            else:
                cnv.append(x)
    
    return cnv

with open(outputFile, 'w') as f:
    print('chromosome', 'start', 'end', 'cn', 'info', 'tool', sep='\t', file=f)
    for x in readFile(cnvkit):
        print(*x, 'cnvkit', sep='\t', file=f)
    for x in readFile(cnvpytor):
        print(*x, 'cnvpytor', sep='\t', file=f)
    for x in readFile(freec):
        print(*x, '-', 'freec', sep='\t', file=f)
    for x in readFile(mops):
        print(*x, '-', 'mops', sep='\t', file=f)