#! /usr/bin/env python
# Merge CNV calling results from 5 tools, assign 'accumulative Score' (1. AS).

import sys
import os

def overlap(cnv1, cnv2):
    '''Calculate overlap proportion of two CNVs, return True if it is larger than some threshold
    - cnv1: ground truth CNV
    - cnv2: detected CNV
    - return: size and maximum proportion of overlapped region between two CNVs
    '''

    c1, s1, e1, cn1 = cnv1[0], int(cnv1[1]), int(cnv1[2]), int(cnv1[3])
    c2, s2, e2, cn2 = cnv2[0], int(cnv2[1]), int(cnv2[2]), int(cnv2[3])
    cnvLen1, cnvLen2 = e1 - s1, e2 - s2
    overlap = 0
    assert cnvLen1 > 0 and cnvLen2 > 0, "The length of CNV is < 0, something wrong!"
    if c1 == c2:
        if (cn1>2 and cn2>2) or (cn1<2 and cn2<2):
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
    
    cnvProp1 = round((overlap/cnvLen1), 2)
    cnvProp2 = round((overlap/cnvLen2), 2)
    cnvProp = max(cnvProp1, cnvProp2)
    return overlap, cnvProp

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
    
    smoove = sys.argv[1]
    delly = sys.argv[2]
    cnvkit = sys.argv[3]
    cnvpytor = sys.argv[4]
    mops = sys.argv[5]
    outputFile = sys.argv[6]

    sample = os.path.basename(smoove).split('.')[0]

    # there is priority for keeping CNVs when merging
    cnvfiles = [smoove, delly, cnvkit, cnvpytor, mops]
    print('Merging CNV results from: ', cnvfiles)
    cnvtools = ['smoove', 'delly', 'cnvkit', 'cnvpytor', 'mops']

    overlapPropThreshold = 0.3
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
        overlapSize = 0
        for i, cnv2 in enumerate(cnvs3):
            # see if two cnvs have overlap
            overlapSize, overlapProp = overlap(cnv1[:4], cnv2[:4])
            if overlapProp == 0:
                continue
            # if yes but overlap proportion is less than a threshold, keep former one and pop out later one
            elif 0 < overlapProp < overlapPropThreshold:
                cnvs2.pop(i-count)
                count += 1
            # if overlap proportion is larger than the threshold, extend breakpoints
            else:
                assert cnv1[-1] != cnv2[-1], "Overlapped CNV from same tool {:s}! Please make sure these conflicts are solved before merging".format(cnv1[-1])
                # as smoove do not determine exact copy number, we use the other tool to supplement
                if tmpCnv[-1] == 'smoove':
                    tmpCnv[3] = cnv2[3]
                tmpCnv[1:3] = [min(cnv1[1], cnv2[1]), max(cnv1[2], cnv2[2])]
                accumLen += overlapSize
                cnvs2.pop(i-count)
                tmpCnv[-1] = ','.join([tmpCnv[-1], cnv2[-1]])
                count += 1
                cnv1 = tmpCnv

        tmpTools = set(cnv1[-1].split(','))
        cnv1.append(len(tmpTools))
        accumFold = round(accumLen * 100 / (int(cnv1[2]) - int(cnv1[1])), 1)   # percentage
        cnv1.append(accumFold)
        
        mergedCnvs.append(cnv1)
        cnvs = cnvs2[:]
        
    with open(outputFile, 'w') as f:
        print('chromosome', 'start', 'end', 'cn', 'tools', 'toolNum', 'accumScore', 'sample', sep='\t', file=f)
        for cnv in mergedCnvs:
            print(*cnv, sample, sep='\t', file=f)

