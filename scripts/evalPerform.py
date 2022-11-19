#! /urs/bin/env python
# Evaluate the performance of different tools by comparing their CNV calling results with simulated
# ground truth set. We use false discovery rate and sensitivity to benchmark.
# Sensitivity = # of true positive CNVs / # of all simulated CNVs
# FDR = # of false positive CNVs / # of all detected CNVs

import sys

def overlap(cnv1, cnv2):
    '''Calculate overlap proportion of two CNVs, return True if it is larger than 0.3'''

    c1, s1, e1, cn1 = cnv1[0], int(cnv1[1]), int(cnv1[2]), int(cnv1[3])
    c2, s2, e2, cn2 = cnv2[0], int(cnv2[1]), int(cnv2[2]), int(cnv2[3])
    cnvLen1 = e1 - s1
    cnvLen2 = e2 - s2
    overlap = 0
    assert cnvLen1 > 0 and cnvLen2 > 0, "The length of CNV is less than 0, wrong!"
    if c1 == c2:
        if (cn1>2 and cn2>2) or (cn1<2 and cn2<2):
            # if s1 < s2 < e1 or s1 < e2 < e1:
            #     print("Truth CNV {:s}:{:d}-{:d} overlaps with called CNV {:s}:{:d}-{:d}".format(c1, 
            #         s1, e1, c2, s2, e2))
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
    if cnvProp1 > 0.3 or cnvProp2 > 0.3:
        return True
    else:
        return False


if __name__ == "__main__":

    callFile = sys.argv[1]      # CNV callset, from merged or single tool
    truthFile = sys.argv[2]     # CNV ground truth file
    Type = sys.argv[3]          # for single-tool or merged result, could be 'single' or 'merge'

    # define some threshold to filter CNV for merged result
    accumScoreThe = 0      # can be adjusted
    goodScoreThe = 90       # can be adjusted
    dupholdScoreThe = 90    # cannot be adjusted

    callCnvs = []
    with open(callFile, 'r') as f:
        for line in f:
            if line.startswith('chrom'):
                continue
            x = line.strip().split('\t')
            cnv = x[:4]
            if Type == 'merge':
                accumScore = float(x[5])
                goodScore = float(x[6])
                dupholdScore = float(x[7])
                if accumScore >= accumScoreThe and goodScore >= goodScoreThe and \
                dupholdScore >= dupholdScoreThe:
                    callCnvs.append(cnv)
            else:
                callCnvs.append(cnv)

    # print('Number of identified CNVs: ', len(callCnvs))

    truthCnvs = []
    with open(truthFile, 'r') as f:
        for line in f:
            if line.startswith('chrom'):
                continue
            x = line.strip().split('\t')
            truthCnvs.append(x)

    # print('Number of truth CNVs: ', len(truthCnvs))

    truthCnvLen = len(truthCnvs)    # number of simulated CNVs
    callCnvLen = len(callCnvs)      # number of called CNVs
    tp = []
    observTP = []
    callCnvs2 = callCnvs[:]
    while truthCnvs:
        cnv1 = truthCnvs.pop(0)
        count = 0       # count the number of CNV in calling set that has overlap with truth CNV
        tmpCnv = []
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

    print("Benchmark for call set {:s} under '{:s}' mode".format(callFile, Type))
    print("Sensitivity is {:.2f}, FDR is {:.2f}".format(sensitivity, fdr))

