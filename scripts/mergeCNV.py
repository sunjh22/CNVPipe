#! /usr/bin/env python
# Merge CNV calling results from 5 tools, assign 'accumulative Score' (1. AS).

import sys
import os
from mergeTruthSet import mergeTruthCNV

def overlapLen(cnvRegion1, cnvRegion2):
    """
    Find the length of overlapped region between two CNVs
    - cnvRegion1: a list contains chromosome, start and end
    - cnvRegion2: a list contains chromosome, start and end
    Return the length of overlapped region
    """

    c1, s1, e1 = cnvRegion1[0], int(cnvRegion1[1]), int(cnvRegion1[2])
    c2, s2, e2 = cnvRegion2[0], int(cnvRegion2[1]), int(cnvRegion2[2])
    cnvLen = e1 - s1
    assert cnvLen > 0, "The length of CNV is less than 0, something wrong! Please check."
    if c1 == c2:
        if s1 < s2 < e1 or s1 < e2 < e1 or s2 < s1 < e1 < e2:
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


def readFile(infile):
    """
    Read list of CNV result from CNV calling tools
    # Return a list of chrom, start, end and cn
    """
    with open(infile, 'r') as f:
        for x in f:
            if x.startswith('chromosome'):
                continue
            x = x.strip().split('\t')
            if not x[0].startswith('chr'):
                x[0] = 'chr' + x[0]
            x = x[:4]
            yield x


if __name__ == "__main__":
    
    cnvkit = sys.argv[1]
    delly = sys.argv[2]
    mops = sys.argv[3]
    cnvpytor = sys.argv[4]
    smoove = sys.argv[5]
    outputFile = sys.argv[6]

    sample = os.path.basename(cnvkit).split('.')[0]

    # there is priority for keeping CNVs when merging
    # cnvfiles = [cnvkit, delly, mops, cnvpytor]
    cnvfiles = [smoove, delly, cnvkit, cnvpytor, mops]
    print(cnvfiles)
    # cnvtools = ['cnvkit', 'delly', 'mops', 'cnvpytor', 'smoove']
    cnvtools = ['smoove', 'delly', 'cnvkit', 'cnvpytor', 'mops']

    cnvs = []
    for i, cnvfile in enumerate(cnvfiles):
        for cnv in readFile(cnvfile):
            cnv.append(cnvtools[i])
            cnvs.append(cnv)    # five columns in cnvs: chr, start, end, cn, tool

    cnvs2 = cnvs[:]
    mergedCnvs = []
    while cnvs:
        cnv1 = cnvs.pop(0)
        tmpCnv = cnvs2.pop(0)
        cnvs3 = cnvs2[:]
        count = 0
        accumLen = 0
        cn1 = int(cnv1[3])
        overlapSize = 0
        for i, cnv2 in enumerate(cnvs3):
            cn2 = int(cnv2[3])
            if (cn1>2 and cn2>2) or (cn1<2 and cn2<2):
                overlapSize = overlapLen(cnv1[:3], cnv2[:3])
                updateCnv1 = mergeTruthCNV(cnv1[:4], cnv2[:4])
                tmpCnv[:4] = updateCnv1
                if overlapSize > 0:
                    assert cnv1[-1] != cnv2[-1], "Overlapped CNV from same tool {:s}! Please make sure these conflicts are solved before merging".format(cnv1[-1])
                    accumLen += overlapSize
                    cnvs2.pop(i-count)
                    tmpCnv[-1] = ','.join([tmpCnv[-1], cnv2[-1]])
                    count += 1
                    cnv1 = tmpCnv

        tmpTools = set(cnv1[-1].split(','))
        cnv1.append(len(tmpTools))
        accumFold = 100 - round(accumLen * 100 / (int(cnv1[2]) - int(cnv1[1])), 1)   # percentage
        cnv1.append(accumFold)
        
        mergedCnvs.append(cnv1)
        cnvs = cnvs2[:]
        
    with open(outputFile, 'w') as f:
        print('chromosome', 'start', 'end', 'cn', 'tools', 'toolNum', 'accumScore', 'sample', sep='\t', file=f)
        for cnv in mergedCnvs:
            print(*cnv, sample, sep='\t', file=f)

