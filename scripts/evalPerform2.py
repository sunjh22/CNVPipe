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

def evaluate(callFile, truthFile, Type):
    '''Calculate sensitivity and FDR for single tool and merged results'''

    # define some threshold to filter CNV for merged result
    accumScoreThe = 0      # can be adjusted
    goodScoreThe = -1000       # can be adjusted
    dupholdScoreThe = 0    # cannot be adjusted
    toolNumThe = 2             # CNVs called by at least how many tools

    callCnvs = []
    with open(callFile, 'r') as f:
        for line in f:
            if line.startswith('chrom'):
                continue
            x = line.strip().split('\t')
            cnv = x[:4]
            if Type == 'merge2':
                accumScore = float(x[5])
                goodScore = float(x[6])
                # dupholdScore = float(x[7])
                toolNum = int(x[4])
                if accumScore >= accumScoreThe and goodScore >= goodScoreThe and \
                toolNum >= toolNumThe:
                    callCnvs.append(cnv)
            else:
                callCnvs.append(cnv)

    truthCnvs = []
    with open(truthFile, 'r') as f:
        for line in f:
            if line.startswith('chrom'):
                continue
            x = line.strip().split('\t')
            truthCnvs.append(x)

    truthCnvLen = len(truthCnvs)    # number of simulated CNVs
    callCnvLen = len(callCnvs)      # number of called CNVs
    tp = []
    observTP = []
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

def evaluateHelper(truthFile, tools, fold, outputFile):
    for tool in tools:
        if tool == 'merge2':
            callFile = '/home/jhsun/data3/project/CNVPipe/analysis/res/' + tool + '/sample' + \
                str(i) + '.bed'
            sensitivity, fdr = evaluate(truthFile=truthFile, callFile=callFile, Type=tool)
            print(fold, 'sample'+str(i), tool, sensitivity, fdr, sep='\t', file=outputFile)
        else:
            callFile = '/home/jhsun/data3/project/CNVPipe/analysis/res/' + tool + '/sample' + \
                str(i) + '.bed'
            sensitivity, fdr = evaluate(truthFile=truthFile, callFile=callFile, Type=tool)
            print(fold, 'sample'+str(i), tool, sensitivity, fdr, sep='\t', file=outputFile)
    print('Finished for {:s} fold.'.format(fold))


if __name__ == "__main__":

    outputFile = sys.argv[1]    # say sensitivity-FDR.v2.txt
    out = open(outputFile, 'w')
    print('fold', 'sample', 'tool', 'sensitivity', 'FDR', sep='\t', file=out)

    tools = ['merge2', 'cnvkit', 'delly', 'mops', 'cnvpytor', 'smoove', 'freec']
    for i in range(1,19):
        truthFile = '/home/jhsun/data3/project/CNVPipe/simulation/simuGenome/tumor' + str(i) + '.truth.bed'
        if 1 <= i <= 6:
            evaluateHelper(truthFile=truthFile, tools=tools, fold='1x', outputFile=out)
        elif 7 <= i <= 12:
            evaluateHelper(truthFile=truthFile, tools=tools, fold='10x', outputFile=out)
        else:
            evaluateHelper(truthFile=truthFile, tools=tools, fold='30x', outputFile=out)
