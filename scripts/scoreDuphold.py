#! /usr/bin/env python
### Score CNV region by DHFFC and DHBFC from duphold

import sys

inputFile = sys.argv[1]
outputFile = sys.argv[2]

with open(inputFile, 'r') as f, open(outputFile, 'w') as g:
    print('chromosome\tstart\tend\tcn\tcnv\taccumScore\tdupholdScore\ttoolName\ttoolNum\tsample', file=g)
    for x in f:
        x = x.strip().split('\t')
        chrom, start, end = x[:3]
        cn = int(x[3])
        aScore = float(x[4])
        dhffc = float(x[5])
        dhbfc = float(x[6])
        toolName = x[7]
        toolNum = x[8]
        spl = x[9]
        if cn < 2:
            cnv = 'deletion'
            if dhffc <= 0.5:
                score = 100
            elif 0.5 < dhffc <= 0.7:
                score = 80
            elif 0.7 < dhffc < 0.99:
                score = 50
            else:
                score = 0
        else:
            cnv = 'duplication'
            if dhbfc >= 1.5:
                score = 100
            elif 1.3 <= dhbfc < 1.5:
                score = 80
            elif 1.01 < dhbfc < 1.3:
                score = 50
            else:
                score = 0
        print(chrom, start, end, cn, cnv, aScore, score, toolName, toolNum, spl, sep='\t', file=g)