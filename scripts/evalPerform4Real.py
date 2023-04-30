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
import pandas as pd
import joblib

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
    # if min(cnvProp1, cnvProp2) > 0.3:
    if cnvProp1 >= 0.5:
        return True
    else:
        return False


def readCNV4Merge(sample_file):
    sample_data_in = pd.read_csv(sample_file, sep='\t', skiprows=1, names=['cn', 'accumScore', 'depthScore', 'tools', 'toolNum', 'cnvfilter', 'goodScore'], usecols=[3,5,6,7,8,9,10])
    sample_data_in = sample_data_in.query('cn < 2')
    sample_data = sample_data_in.copy()
    sample_data.loc[:, 'cnvfilter'] = [0 if x==True else 1 for x in sample_data_in['cnvfilter']]
    toolScore = {'smoove': 5, 'delly': 4, 'cnvkit': 3, 'cnvpytor': 2, 'mops': 1}
    sample_data.loc[:, 'tools'] = [sum(toolScore[t] for t in set(x.split(','))) for x in sample_data_in['tools']]
    sample_data = sample_data.drop(columns='cn')
    return sample_data


def evaluate(truthFile, callFile, Type, deletion=True):
    '''Calculate sensitivity and FDR for single tool and merged results. We generally take CNVs
    detected by more than one tool in merged results as input.
    - truthFile: a file with ground truth CNVs
    - callFile: a file with detected CNVs
    - Type: 'merge' or other tools
    - dup: whether the truth set of sample has duplications
    - return: sensitivity and FDR
    '''

    # read detected CNVs, read more columns for CNV filtering for merged CNV set
    callCnvs = []
    with open(callFile, 'r') as f:
        for line in f:
            if line.startswith('chrom'):
                continue
            x = line.strip().split('\t')
            cnv = x[:4]
            cn = int(cnv[3])
            if deletion and cn > 2:
                continue
            # callCnvs.append(cnv)

            # brute-force method
            # if (toolNum >= 2 or 'smoove' in toolName or 'delly' in toolName) and dupholdScore >=50 and cnvfilter == 'True' and goodScore == 100:
            # if 'smoove' in toolName and dupholdScore>0:
            # if ('smoove' in toolName or 'delly' in toolName) and dupholdScore > 0:
            # if toolNum >= 2 and dupholdScore >= 90 and cnvfilter == 'True' and goodScore == 100:
            if Type == 'CNVPipe':
                dupholdScore = int(x[6])
                toolName = x[7].split(',')
                toolNum = int(x[8])
                cnvfilter = x[9]
                goodScore = int(x[10])
                if ('smoove' in toolName or 'delly' in toolName) and cnvfilter == 'True' and goodScore == 100:
                    callCnvs.append(cnv)
            else:
                callCnvs.append(cnv)

    # SVM method
    # if Type == 'merge':
    #     # print("Number of copy number deletions:", len(callCnvs))
    #     tmp_callCnvs = []
    #     clf = joblib.load("/data3/jhsun/github-repo/CNVPipe/resources/SVM/cnv_svm_classifier_real_10x.pkl")
    #     allCnvs = readCNV4Merge(sample_file=callFile)
    #     predictions = clf.predict(allCnvs)
    #     # print("Number of CNVs input into SVM:", len(predictions))
    #     for idx, label in enumerate(predictions):
    #         if label == 'T':
    #             tmp_callCnvs.append(callCnvs[idx])
    #     callCnvs = tmp_callCnvs[:]
    
    # read ground truth CNVs
    truthCnvs = []
    with open(truthFile, 'r') as f:
        for line in f:
            if line.startswith('chrom'):
                continue
            x = line.strip().split('\t')[:4]
            if deletion:
                if x[3] not in ['deletion', 'DEL']:
                    continue
                x[3] = 1
            else:
                if x[3] in ['deletion', 'DEL']:
                    x[3] = 1
                elif x[3] in ['duplication', 'DUP']:
                    x[3] = 3
                else:
                    pass
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

    sensitivity = round(len(tp)/truthCnvLen, 3)
    if callCnvLen == 0:
        return sensitivity, 1
    fdr = round((callCnvLen-len(observTP))/callCnvLen, 3)
    return sensitivity, fdr


def evaluateHelper(truthFile, tools, sampleID, outputFile):
    for tool in tools:
        callFile = '/home/jhsun/data3/project/CNVPipe/realAnalysis/res/' + tool + '/' + \
                sampleID + '.bed'
        if sampleID in ['sample13']:
            sensitivity, fdr = evaluate(truthFile=truthFile, callFile=callFile, Type=tool, deletion=False)
        else:
            callFile = '/home/jhsun/data3/project/CNVPipe/realAnalysis/res/' + tool + '/' + \
                sampleID + '.bed'
            sensitivity, fdr = evaluate(truthFile=truthFile, callFile=callFile, Type=tool, deletion=True)
        
        print(sampleID, tool, sensitivity, fdr, sep='\t')
        print(sampleID, tool, sensitivity, fdr, sep='\t', file=outputFile)



if __name__ == "__main__":

    outputFile = sys.argv[1]
    out = open(outputFile, 'w')
    print('sample', 'tool', 'sensitivity', 'FDR', sep='\t')
    print('sample', 'tool', 'sensitivity', 'FDR', sep='\t', file=out)

    tools = ['CNVPipe', 'cnvkit', 'delly', 'cnvpytor', 'smoove', 'mops']
    samples = ['NA12878-1', 'NA12878-2', 'CHM13', 'AK1', 'HG002', 'HG00514', 'HG00733', 'NA19240']

    for sampleID in samples:
        truthFile = '/home/jhsun/data3/project/CNVPipe/realAnalysis/truthSet/' + sampleID + '-SVset.1kb.bed'
        evaluateHelper(truthFile=truthFile, tools=tools, sampleID=sampleID, outputFile=out)
