#! /usr/bin/env python
# This merge script is used for median-depth data (between 0.5x and 5x depth).
# We will merge CNV calling results from mops, cnvkit, delly, cnvpytor and smoove. 
# Assign 'accumulative Score'(1. AS).

import sys
import os
from mergeCNV import overlap

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
    
    mops = sys.argv[1]
    cnvkit = sys.argv[2]
    delly = sys.argv[3]
    cnvpytor = sys.argv[4]
    smoove = sys.argv[5]
    outputFile = sys.argv[6]

    sample = os.path.basename(mops).split('.')[0]

    # there is priority for keeping CNVs when merging
    cnvfiles = [mops, cnvkit, delly, cnvpytor, smoove]
    print('Merging CNV results from: ', cnvfiles)
    cnvtools = ['mops', 'cnvkit', 'delly', 'cnvpytor', 'smoove']

    overlapPropThreshold = 0.3
    cnvs = []
    for i, cnvfile in enumerate(cnvfiles):
        for cnv in readFile(cnvfile):
            cnv.append(cnvtools[i])
            cnvs.append(cnv)    # five columns in cnvs: chr, start, end, cn, tool

    cnvs2 = cnvs[:]
    mergedCnvs = []
    # we designed a strategy to recursively merging overlapped CNVs and extending their breakpoints
    while cnvs:
        cnv1 = cnvs.pop(0)
        tmpCnv = cnvs2.pop(0)
        cnvs3 = cnvs2[:]
        count = 0
        accumLen = 0
        cn1 = int(cnv1[3])
        overlapSize = 0
        for i, cnv2 in enumerate(cnvs3):
            # see if two cnvs have overlap
            overlapSize, overlapProp = overlap(cnv1[:4], cnv2[:4])
            if overlapProp == 0:
                continue
            # if yes but overlap proportion is less than 0.5, keep former one and pop out later one
            elif 0 < overlapProp < overlapPropThreshold:
                cnvs2.pop(i-count)
                count += 1
            # if overlap proportion is larger than 0.5, extend breakpoints
            else:
                assert cnv1[-1] != cnv2[-1], "Overlapped CNV from same tool {:s}! Please make sure these conflicts are solved before merging".format(cnv1[-1])
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

