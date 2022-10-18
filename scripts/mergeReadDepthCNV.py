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


def overlapBad(region, bad):
    # filter CNV region if 30% of it overlaps with known low-map region
    c1, s1, e1 = region[0], int(region[1]), int(region[2])
    c2, s2, e2 = bad[0], int(bad[1]), int(bad[2])
    regionSize = e1 - s1
    badSize = e2 - s2
    assert regionSize > 0, "The length of CNV region is 0, something wrong happens!"
    assert badSize > 0, "The length of bad region is 0, something wrong happens!"
    if c1 == c2:
        if s2 <= s1 <= e2:
            overlap = e2 - s1
            if overlap > regionSize * 0.3:
                return True
            else:
                return False
        if s2 <= e1 <= e2:
            overlap = e1 - s2
            if overlap > regionSize * 0.3:
                return True
            else:
                return False


def readFile(infile, bad):
    with open(infile, 'r') as f:
        for x in f:
            flag = False
            if x.startswith('chromosome'):
                continue
            x = x.strip().split('\t')
            if not x[0].startswith('chr'):
                x[0] = 'chr' + x[0]

            # make cnv results from different tools output in the same format
            if len(x) > 5:
                x.pop(4)
            elif len(x) == 5:
                x.pop(4)
                x.append('-')
            else:
                x.append('-')
            
            # prepare for CNVfilteR input
            if int(x[3]) > 2:
                x.append('duplication')
            else:
                x.append('deletion')

            # remove CNVs in low-map regions
            for y in bad:
                if overlapBad(x[:3], y):
                    flag = True
                    break
            if flag:
                break
            else:
                yield x


if __name__ == "__main__":
    
    cnvkit = sys.argv[1]
    cnvpytor = sys.argv[2]
    freec = sys.argv[3]
    mops = sys.argv[4]
    lowMapFile = sys.argv[5]    # by default, we use blacklist from 10x
    outputFile = sys.argv[6]

    sample = os.path.basename(cnvkit).split('.')[0]

    bad = readBadRegion(lowMapFile)

    with open(outputFile, 'w') as f:
        print('chromosome', 'start', 'end', 'cn', 'info', 'cnv', 'sample', 'tool', sep='\t', file=f)
        for x in readFile(cnvkit, bad):
            print(*x, sample, 'cnvkit', sep='\t', file=f)
        for x in readFile(cnvpytor, bad):
            print(*x, sample, 'cnvpytor', sep='\t', file=f)
        for x in readFile(freec, bad):
            print(*x, sample, 'freec', sep='\t', file=f)
        for x in readFile(mops, bad):
            print(*x, sample, 'mops', sep='\t', file=f)