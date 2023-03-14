#! /usr/bin/env python

import sys

def overlap(cnv1, cnv2):
    '''Calculate overlap proportion of two CNVs, return True if it is larger than some threshold
    - cnv1: ground truth CNV
    - cnv2: detected CNV
    - return: size of overlapped region between two CNVs
    '''

    c1, s1, e1 = cnv1[0], int(cnv1[1]), int(cnv1[2])
    c2, s2, e2 = cnv2[0], int(cnv2[1]), int(cnv2[2])
    cnvLen1, cnvLen2 = e1 - s1, e2 - s2
    overlap = 0
    assert cnvLen1 > 0 and cnvLen2 > 0, "The length of CNV is < 0, something wrong!"
    if c1 == c2:
        if s2 <= s1 < e1 <= e2:
            overlap = e1 - s1
        elif s2 <= s1 <= e2 < e1:
            overlap = e2 - s1
        elif s1 <= s2 < e2 <= e1:
            overlap = e2 - s2
        elif s1 < s2 <= e1 <= e2:
            overlap = e1 - s2
        else:
            pass

    return overlap

def readRegionFile(infile):
    bad = []
    with open(infile, 'r') as f:
        for x in f:
            x = x.strip().split('\t')[:3]
            bad.append(x)
    return bad

def calculateOverlapScore(target_cnv, cnv_list, reciprocal_prop):
    """
    Calculate the accumulative overlap between a target CNV and a list of cnvs, only cnvs (in cnv 
    list) with over 30% reciprocal overlap with target CNV are selected. Reduce the score 100 by the
    overlap fraction and the number of overlapped cnvs.
    - target_cnv: a single CNV with chromosome, start and end
    - cnv_list: a list of cnv regions. Could be a list of bad regions or common normal CNVs
    Return a score, with 100 means no overlap with cnv list at all.
    """
    accumLen = 0    # accumulative length of overlap region between cnv and bad
    i = 0
    print("A new CNV region: ", target_cnv)

    targetSize = int(target_cnv[2]) - int(target_cnv[1])

    for cnv in cnv_list:
        cnvSize = int(cnv[2]) - int(cnv[1])
        overlapSize = overlap(target_cnv, cnv)
        overlapProp = max(overlapSize/targetSize, overlapSize/cnvSize)
        if overlapProp > reciprocal_prop:
            i += 1
            accumLen += overlapSize

    accumProp = round(accumLen * 100 / (int(target_cnv[2]) - int(target_cnv[1])))
    score = 100 - accumProp - i * 10
    print("This CNV totally overlaps with {:d} CNV regions\n".format(i))
    return score


if __name__ == "__main__":
    
    inputFile = sys.argv[1]
    lowMapFile = sys.argv[2]    # by default, we use blacklist from 10x
    normalCNVFile = sys.argv[3] # common SV list in normal population
    outputFile = sys.argv[4]

    bad_list = readRegionFile(lowMapFile)
    normal_list = readRegionFile(normalCNVFile)

    with open(inputFile, 'r') as f, open(outputFile, 'w') as g:
        for x in f:
            if x.startswith('chromosome'):
                continue
            cnv = x.strip().split('\t')[:3]
            bad_score = calculateOverlapScore(cnv, bad_list, 0.3)
            normal_score = calculateOverlapScore(cnv, normal_list, 0.5)
            print(x.strip(), bad_score, normal_score, sep='\t', file=g)
    