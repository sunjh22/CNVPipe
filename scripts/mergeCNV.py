#! /usr/bin/env python
# Merge CNV calling results from six tools, assign accumScore.
# And calculate overlap proportion with SV black regions, assign goodScore

import sys
import os

def readBadRegion(infile):
    bad = []
    with open(infile, 'r') as f:
        for x in f:
            x = x.strip().split('\t')[:3]
            bad.append(x)
    return bad

def overlapLen(cnvRegion, badRegion):
    c1, s1, e1 = cnvRegion[0], int(cnvRegion[1]), int(cnvRegion[2])
    c2, s2, e2 = badRegion[0], int(badRegion[1]), int(badRegion[2])
    cnvLen = e1 - s1
    assert cnvLen > 0, "The length of CNV is less than 0, wrong!"
    if c1 == c2:
        if s1 < s2 < e1 or s1 < e2 < e1:
            print("CNV region {:s}:{:d}-{:d} overlaps with bad region {:s}:{:d}-{:d}".format(c1, 
                s1, e1, c2, s2, e2))
        if s2 <= s1 < e1 <= e2:
            overlap = e1 - s1
            return overlap
        elif s2 <= s1 <= e2 < e1:
            overlap = e2 - s1
            return overlap
        elif s1 <= s2 < e2 <= e1:
            overlap = e2 - s2
            return overlap
        elif s1 < s2 <= e1 <= e2:
            overlap = e1 - s2
            return overlap
        else:
            return 0
    
    return 0

def overlapScore(cnvRegion, bad):
    accumLen = 0    # accumulative length of overlap region between cnv and bad
    i = 0
    print("A new CNV region: ", cnvRegion)

    for badRegion in bad:
        overlapSize = overlapLen(cnvRegion, badRegion)
        if overlapSize > 0:
            i += 1
            accumLen += overlapSize

    accumProp = round(accumLen * 100 / (int(cnvRegion[2]) - int(cnvRegion[1])))
    score = 100 - accumProp - i * 10
    print("This CNV totally overlaps with {:d} bad genomic regions\n".format(i))
    return score

def readFile(infile):
    with open(infile, 'r') as f:
        for x in f:
            if x.startswith('chromosome'):
                continue
            x = x.strip().split('\t')
            if not x[0].startswith('chr'):
                x[0] = 'chr' + x[0]
            # only reserve chrom, start, end and cn
            x = x[:4]
            yield x


if __name__ == "__main__":
    
    cnvkit = sys.argv[1]
    delly = sys.argv[2]
    smoove = sys.argv[3]
    mops = sys.argv[4]
    cnvpytor = sys.argv[5]
    # freec = sys.argv[6]
    lowMapFile = sys.argv[6]    # by default, we use blacklist from 10x
    outputFile = sys.argv[7]

    sample = os.path.basename(cnvkit).split('.')[0]

    # there is priority for keeping CNVs when merging
    cnvfiles = [cnvkit, delly, smoove, mops, cnvpytor]
    print(cnvfiles)
    cnvtools = ['cnvkit', 'delly', 'smoove', 'mops', 'cnvpytor']
    bad = readBadRegion(lowMapFile)

    cnvs = []
    for i, cnvfile in enumerate(cnvfiles):
        for cnv in readFile(cnvfile):
            cnv.append(cnvtools[i])
            cnvs.append(cnv)    # five columns in cnvs: chr, start, end, cn, tool

    tmpCnvs = cnvs[:]
    mergedCnvs = []
    while cnvs:
        cnv1 = cnvs.pop(0)
        tmpCnv = tmpCnvs.pop(0)
        accumLen = 0
        cn1 = int(cnv1[3])
        count = 0   # count the number of cnv2 that overlaps with cnv1
        for i, cnv2 in enumerate(cnvs):
            cn2 = int(cnv2[3])
            if (cn1>2 and cn2>2) or (cn1<2 and cn2<2):
                overlapSize = overlapLen(cnv1[:3], cnv2[:3])
                if overlapSize > 0:
                    assert cnv1[-1] != cnv2[-1], "Overlapped CNV from same tool {:s}!".format(cnv1[-1])
                    accumLen += overlapSize
                    tmpCnvs.pop(i-count)
                    tmpCnv[-1] = ','.join([tmpCnv[-1], cnv2[-1]])
                    count += 1

        tmpTools = set(tmpCnv[-1].split(','))
        tmpCnv[-1] = len(tmpTools)
        accumFold = round(accumLen * 100 / (int(cnv1[2]) - int(cnv1[1])), 1)   # percentage
        tmpCnv.append(accumFold)
        
        mergedCnvs.append(tmpCnv)   # six columns: chr, start, end, cn, tools, accumFold
        cnvs = tmpCnvs[:]
        
    with open(outputFile, 'w') as f:
        print('chromosome', 'start', 'end', 'cn', 'toolNum', 'accumScore', 'goodScore', 
            'sample', sep='\t', file=f)
        for cnv in mergedCnvs:
            # score CNVs by calculating the overlap proportion with bad regions
            score = overlapScore(cnv[:3], bad)
            print(*cnv, score, sample, sep='\t', file=f)

