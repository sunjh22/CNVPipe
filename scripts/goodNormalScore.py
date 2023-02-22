#! /usr/bin/env python

import sys
from mergeCNV3 import overlapLen

def readBadRegion(infile):
    bad = []
    with open(infile, 'r') as f:
        for x in f:
            x = x.strip().split('\t')[:3]
            bad.append(x)
    return bad

def overlapScore4BadRegion(cnvRegion, bad_list):
    """
    Calculate the accumulative overlap between one CNV and a list of bad genomic regions, use a
    score to represent. Both the overlap proportion and the number of overlaps count.
    - cnvRegion: a single CNV with chromosome, start and end
    - bad: a list of genomic bad regions
    Return a score with 100 means no overlap with bad region at all.
    """
    accumLen = 0    # accumulative length of overlap region between cnv and bad
    i = 0
    print("A new CNV region: ", cnvRegion)

    for badRegion in bad_list:
        overlapSize = overlapLen(cnvRegion, badRegion)
        if overlapSize > 0:
            i += 1
            accumLen += overlapSize

    accumProp = round(accumLen * 100 / (int(cnvRegion[2]) - int(cnvRegion[1])))
    score = 100 - accumProp - i * 10
    print("This CNV totally overlaps with {:d} bad genomic regions\n".format(i))
    return score

def overlapScore4NormalPopu(cnvRegion, normal_list):
    pass

if __name__ == "__main__":
    
    inputFile = sys.argv[1]
    lowMapFile = sys.argv[2]    # by default, we use blacklist from 10x
    normalCNVFile = sys.argv[3]
    outputFile = sys.argv[4]

    bad_list = readBadRegion(lowMapFile)
    normal_list = []

    with open(inputFile, 'r') as f, open(outputFile, 'w') as g:
        for x in f:
            cnv = x.strip().split('\t')[:3]
            score1 = overlapScore4BadRegion(cnv, bad_list)
            score2 = overlapScore4NormalPopu(cnv, normal_list)
            print(x, score1, score2, sep='\t', file=g)
    