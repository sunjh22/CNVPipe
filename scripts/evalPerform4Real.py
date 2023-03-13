#! /urs/bin/env python
# Usage: python evalPerform4Real.py 

# Evaluate the performance of different tools for real WGS data by comparing their CNV calling 
# results with ground truth set. We use false discovery rate and sensitivity to benchmark.

# Generally, we define true positive calls as CNV has 50% reciprocal overlap with truth set. For
# CNVPipe results, we only take those CNVs detected by more than 
# one tool as the input.

# Sensitivity = # of true positive CNVs / # of all simulated CNVs
# FDR = # of false positive CNVs / # of all detected CNVs
# In this case, recall = sensitivity, precision = 1 - FDR.

import sys
import os

def overlap(cnv1, cnv2):
    '''Calculate overlap proportion of two CNVs, return True if it is larger than some threshold
    - cnv1: ground truth CNV
    - cnv2: detected CNV
    - return: True or False for whether two CNVs overlapped
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
    # if cnvProp1 > 0.8 and cnvProp2 > 0.8:
    if cnvProp1 >= 0.5:
        return True
    else:
        return False


def evaluate(truthFile, callFile, Type):
    '''Calculate sensitivity and FDR for single tool and merged results. We generally take CNVs
    detected by more than one tool in merged results as input.
    - truthFile: a file with ground truth CNVs
    - callFile: a file with detected CNVs
    - return: sensitivity and FDR
    '''

    # define some threshold to filter CNV for merged result
    dupholdScoreThe = 90    # cannot be adjusted
    toolNumThe = 2         # CNVs called by at least how many tools, equal to accumScoreThe > 0

    # read detected CNVs, read more columns for CNV filtering for merged CNV set
    callCnvs = []
    with open(callFile, 'r') as f:
        for line in f:
            if line.startswith('chrom'):
                continue
            x = line.strip().split('\t')
            cnv = x[:4]
            if Type == 'merge':
                dupholdScore = float(x[6])
                toolName = x[7].split(',')
                toolNum = int(x[8])
                cnvfilter = x[9]
                if (toolNum >= toolNumThe or 'smoove' in toolName) and dupholdScore >= dupholdScoreThe and cnvfilter == 'True':
                # if toolNum > 2 or (toolNum==2 and dupholdScore > dupholdScoreThe):
                    callCnvs.append(cnv)
            else:
                callCnvs.append(cnv)

    # read ground truth CNVs
    truthCnvs = []
    with open(truthFile, 'r') as f:
        for line in f:
            if line.startswith('chrom'):
                continue
            x = line.strip().split('\t')[:4]
            if x[3] == 'deletion' or x[3] == 'DEL':
                x[3] = 1
            else:
                x[3] = 3
            truthCnvs.append(x)

    truthCnvLen = len(truthCnvs)    # number of simulated CNVs
    callCnvLen = len(callCnvs)      # number of called CNVs

    # the reason for tp and observTP is that several called CNVs might overlap with the same truth CNV.
    tp = []     # CNVs in truth set that can be identified by calling
    observTP = []       # CNVs in called set that have overlap with truth-set CNVs
    callCnvs2 = callCnvs[:]
    while truthCnvs:
        cnv1 = truthCnvs.pop(0)
        count = 0       # count the number of CNV in calling set that has overlap with truth CNV
        flag = 0        # mark whether the CNV in truth set has been called
        for i, cnv2 in enumerate(callCnvs):
            if overlap(cnv1, cnv2):
                observTP.append(callCnvs2.pop(i-count))
                count += 1
                flag = 1

        if flag == 1:
            tp.append(cnv1)
        callCnvs = callCnvs2[:]

    sensitivity = round(len(tp)/truthCnvLen, 2)
    fdr = round((callCnvLen-len(observTP))/callCnvLen, 2)
    return(sensitivity, fdr)


def evaluateHelper(truthFile, tools, sampleID, outputFile):
    for tool in tools:
        if tool == 'merge':
            callFile = '/home/jhsun/data3/project/CNVPipe/realAnalysis/res/' + tool + '/' + \
                sampleID + '.final.bed'
            sensitivity, fdr = evaluate(truthFile=truthFile, callFile=callFile, Type=tool)
            print(sampleID, tool, sensitivity, fdr, sep='\t')
            print(sampleID, tool, sensitivity, fdr, sep='\t', file=outputFile)
        else:
            callFile = '/home/jhsun/data3/project/CNVPipe/realAnalysis/res/' + tool + '/' + \
                sampleID + '.bed'
            sensitivity, fdr = evaluate(truthFile=truthFile, callFile=callFile, Type=tool)
            print(sampleID, tool, sensitivity, fdr, sep='\t')
            print(sampleID, tool, sensitivity, fdr, sep='\t', file=outputFile)


if __name__ == "__main__":

    outputFile = sys.argv[1]
    out = open(outputFile, 'w')
    print('sample', 'tool', 'sensitivity', 'FDR', sep='\t')
    print('sample', 'tool', 'sensitivity', 'FDR', sep='\t', file=out)

    tools = ['merge', 'cnvkit', 'delly', 'cnvpytor', 'smoove', 'mops']
    samples = ['NA12878-1', 'NA12878-2', 'CHM13', 'AK1', 'HG002', 'HG00514', 'HG00733', 'NA19240', 'sample13', 'sample14']

    for sampleID in samples:
        truthFile = '/home/jhsun/data3/project/CNVPipe/realAnalysis/truthSet/' + sampleID + '-SVset.1kb.bed'
        evaluateHelper(truthFile=truthFile, tools=tools, sampleID=sampleID, outputFile=out)
