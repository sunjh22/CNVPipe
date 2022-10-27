#! /usr/bin/env python

import sys
import os

def readBadRegion(infile):
    bad = []
    with open(infile, 'r') as f:
        for x in f:
            x = x.strip().split('\t')[:3]
            bad.append(x)
    return bad

def overlapProp(cnvRegion, badRegion):
    c1, s1, e1 = cnvRegion[0], int(cnvRegion[1]), int(cnvRegion[2])
    c2, s2, e2 = badRegion[0], int(badRegion[1]), int(badRegion[2])
    cnvLen = e1 - s1
    assert cnvLen > 0, "The length of CNV is less than 0, wrong!"
    if c1 == c2:
        if s2 <= s1 <= e2:
            print("CNV region {:s}:{:d}-{:d} overlaps with bad region {:s}:{:d}-{:d}".format(c1, 
                s1, e1, c2, s2, e2))
            overlap = round((e2 - s1) * 100 / cnvLen)
            return overlap
        if s2 <= e1 <= e2:
            print("CNV region {:s}:{:d}-{:d} overlaps with bad region {:s}:{:d}-{:d}".format(c1, 
                s1, e1, c2, s2, e2))
            overlap = round((e1 - s2) * 100 / cnvLen)
            return overlap
    
    return 0

def overlapScore(cnvRegion, bad):
    accumLen = 0    # accumulative length of overlap region between cnv and bad
    i = 0
    print("A new CNV region: ", cnvRegion)

    for badRegion in bad:
        overlapSize = overlapProp(cnvRegion, badRegion)
        if overlapSize > 0:
            i += 1
            accumLen += overlapSize

    score = 100 - accumLen
    print("This CNV totally overlaps with {:d} bad genomic regions\n".format(i))
    return score

def readFile(infile, bad):
    with open(infile, 'r') as f:
        for x in f:
            if x.startswith('chromosome'):
                continue
            x = x.strip().split('\t')
            if not x[0].startswith('chr'):
                x[0] = 'chr' + x[0]

            # only reserve chrom, start, end and cn
            x = x[:4]

            # score CNVs by calculating the overlap proportion with bad regions
            score = overlapScore(x[:3], bad)
            x.append(str(score))
            yield x


if __name__ == "__main__":
    
    cnvkit = sys.argv[1]
    cnvpytor = sys.argv[2]
    freec = sys.argv[3]
    mops = sys.argv[4]
    smoove = sys.argv[5]
    delly = sys.argv[6]
    lowMapFile = sys.argv[7]    # by default, we use blacklist from 10x
    outputFile = sys.argv[8]

    sample = os.path.basename(cnvkit).split('.')[0]

    cnvfiles = [cnvkit, cnvpytor, freec, mops, smoove, delly]
    cnvtools = ['cnvkit', 'cnvpytor', 'freec', 'mops', 'smoove', 'delly']
    bad = readBadRegion(lowMapFile)

    with open(outputFile, 'w') as f:
        print('chromosome', 'start', 'end', 'cn', 'badRegionScore', 'sample', 'tool', sep='\t', file=f)
        for i, cnvfile in enumerate(cnvfiles):
            for cnv in readFile(cnvfile, bad):
                print(*cnv, sample, cnvtools[i], sep='\t', file=f)