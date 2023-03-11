#! /usr/bin/env python

import sys
from collections import defaultdict

normalBedFile = sys.argv[1]
pathoBedFile = sys.argv[2]
outputFile = sys.argv[3]

pathoCNV = defaultdict(list)
with open(pathoBedFile, 'r') as f:
    for x in f:
        x = x.strip().split('\t')
        cnv = "_".join(x[1:4])
        pathoInfo = [x[5], x[6], x[-2]]
        pathoCNV[cnv] = pathoInfo

with open(normalBedFile, 'r') as f, open(outputFile, 'w') as g:
    print('chrom\tstart\tend\tcn\ttype\tAS\tDS\ttoolName\ttoolNum\tcnvfilter\tGS\tNS\tpathogenicity\tpathoScore\tdosageGene', file=g)
    for x in f:
        if x.startswith('chromosome'):
            continue
        fileds = x.strip().split('\t')
        cnv = "_".join(fileds[:3])
        if cnv in pathoCNV.keys():
            print(x.strip(), *pathoCNV[cnv], sep='\t', file=g)
        else:
            print("The pathogenicity of CNV {:s} was not predicted, please check".format(cnv))
