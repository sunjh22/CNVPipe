#! /urs/bin/env python
# Usage: python evalPerform v1.txt

# Evaluate the performance of different tools by comparing their CNV calling results with simulated
# ground truth set. We use false discovery rate and sensitivity to benchmark.
# Generally, we define true positive calls as CNV has 50% or 80% reciprocal overlap with truth. For
# CNVPipe results of high-depth data (10x and 30x), we only take those CNVs detected by more than 
# one tool as the input. While for CNVPipe results of low-depth data (0.1x and 0.5x), we will 
# take cn.mops results as priority as from our simulation result cn.mops performs better than other
# tools for low-depth data.

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
    # if min(cnvProp1, cnvProp2) > 0.5:
    if cnvProp1 > 0.5:
        return True
    else:
        return False


def readCNV4Merge(sample_file):
    sample_data_in = pd.read_csv(sample_file, sep='\t', skiprows=1, names=['accumScore', 'depthScore', 'tools', 'toolNum', 'cnvfilter'], usecols=[5,6,7,8,9])
    sample_data = sample_data_in.copy()
    sample_data.loc[:, 'cnvfilter'] = [0 if x==True else 1 for x in sample_data_in['cnvfilter']]
    toolScore = {'smoove': 5, 'delly': 4, 'cnvkit': 3, 'cnvpytor': 2, 'mops': 1}
    sample_data.loc[:, 'tools'] = [sum(toolScore[t] for t in set(x.split(','))) for x in sample_data_in['tools']]
    return sample_data


def evaluate(truthFile, callFile, Type, fold):
    '''Calculate sensitivity and FDR for single tool and merged results. We generally take CNVs
    detected by more than one tool in merged results as input.
    - truthFile: a file with ground truth CNVs
    - callFile: a file with detected CNVs
    - return: sensitivity and FDR
    '''

    # define some threshold to filter CNV for merged result
    # accumScoreThe = 0      # can be adjusted
    dupholdScoreThe = 30    # cannot be adjusted
    toolNumThe = 2         # CNVs called by at least how many tools, equal to accumScoreThe > 0

    # read detected CNVs, read more columns for CNV filtering for merged CNV set
    callCnvs = []
    with open(callFile, 'r') as f:
        for line in f:
            if line.startswith('chrom'):
                continue
            x = line.strip().split('\t')
            cnv = x[:4]
            callCnvs.append(cnv)
            # if Type == 'merge':
            #     # accumScore = float(x[5])
            #     # dupholdScore = float(x[6])
            #     # toolNum = int(x[7])
            #     toolNum = int(x[5])
            #     if toolNum >= toolNumThe:
            #     # if toolNum > 2 or (toolNum==2 and dupholdScore > dupholdScoreThe):
            #         callCnvs.append(cnv)
            # else:
            #     callCnvs.append(cnv)

    # SVM method
    if Type == 'merge':
        print("Number of copy number deletions:", len(callCnvs))
        tmp_callCnvs = []
        clf = joblib.load("/data3/jhsun/github-repo/CNVPipe/resources/SVM/cnv_svm_classifier_simu_"+fold+".pkl")
        allCnvs = readCNV4Merge(sample_file=callFile)
        predictions = clf.predict(allCnvs)
        print("Number of CNVs input into SVM:", len(predictions))
        for idx, label in enumerate(predictions):
            if label == 'T':
                tmp_callCnvs.append(callCnvs[idx])
        callCnvs = tmp_callCnvs[:]

    # read ground truth CNVs
    truthCnvs = []
    with open(truthFile, 'r') as f:
        for line in f:
            if line.startswith('chrom'):
                continue
            x = line.strip().split('\t')
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


def evaluate4LowDepth(truthFile, callFile, Type, fold):
    '''Calculate sensitivity and FDR for single tool and merged results. For low read-depth result,
    we take cn.mops result as first priority for merged results
    - truthFile: a file with ground truth CNVs
    - callFile: a file with detected CNVs
    - return: sensitivity and FDR
    '''

    # define some threshold to filter CNV for merged result
    # accumScoreThe = 0      # can be adjusted
    dupholdScoreThe = 30    # cannot be adjusted
    toolNumThe = 2         # CNVs called by at least how many tools, equal to accumScoreThe > 0

    # read detected CNVs, read more columns for CNV filtering for merged CNV set
    callCnvs = []
    with open(callFile, 'r') as f:
        for line in f:
            if line.startswith('chrom'):
                continue
            x = line.strip().split('\t')
            cnv = x[:4]
            callCnvs.append(cnv)
            # if Type == 'merge':

            #     # for sample.final.bed
            #     # accumScore = float(x[5])
            #     # dupholdScore = float(x[6])
            #     # toolNum = int(x[7])

            #     # for sample.merged.bed
            #     toolNum = int(x[5])
            #     toolName = x[4].split(',')
            #     if toolNum >= toolNumThe or 'mops' in toolName:
            #     # if toolNum > 2 or (toolNum==2 and dupholdScore > dupholdScoreThe):
            #         callCnvs.append(cnv)
            # else:
            #     callCnvs.append(cnv)
    
    # SVM method
    if Type == 'merge':
        print("Number of copy number deletions:", len(callCnvs))
        tmp_callCnvs = []
        clf = joblib.load("/data3/jhsun/github-repo/CNVPipe/resources/SVM/cnv_svm_classifier_simu_"+fold+".pkl")
        allCnvs = readCNV4Merge(sample_file=callFile)
        predictions = clf.predict(allCnvs)
        print("Number of CNVs input into SVM:", len(predictions))
        for idx, label in enumerate(predictions):
            if label == 'T':
                tmp_callCnvs.append(callCnvs[idx])
        callCnvs = tmp_callCnvs[:]

    # read ground truth CNVs
    truthCnvs = []
    with open(truthFile, 'r') as f:
        for line in f:
            if line.startswith('chrom'):
                continue
            x = line.strip().split('\t')
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


def evaluateHelper(truthFile, tools, fold, outputFile, sampleID):
    for tool in tools:
        if tool == 'merge':
            callFile = '/home/jhsun/data3/project/CNVPipe/analysis-CNVSimulator/res/' + tool + '/sample' + \
                sampleID + '.bed'
            if fold in ['0.1x', '0.5x']:
                sensitivity, fdr = evaluate4LowDepth(truthFile=truthFile, callFile=callFile, Type=tool, fold=fold)
            else:
                sensitivity, fdr = evaluate(truthFile=truthFile, callFile=callFile, Type=tool, fold=fold)
            print(fold, 'sample'+sampleID, tool, sensitivity, fdr, sep='\t', file=outputFile)
        else:
            callFile = '/home/jhsun/data3/project/CNVPipe/analysis-CNVSimulator/res/' + tool + '/sample' + \
                sampleID + '.bed'
            if not os.path.exists(callFile):
                print(fold, 'sample'+sampleID, tool, 0.001, 0.001, sep='\t', file=outputFile)
                continue
            sensitivity, fdr = evaluate(truthFile=truthFile, callFile=callFile, Type=tool, fold=fold)
            print(fold, 'sample'+sampleID, tool, sensitivity, fdr, sep='\t', file=outputFile)
    print('Finished for {:s} fold.'.format(fold))


if __name__ == "__main__":

    outputFile = sys.argv[1]    # say v1.txt
    out = open(outputFile, 'w')
    print('fold', 'sample', 'tool', 'sensitivity', 'FDR', sep='\t', file=out)

    tools = ['merge', 'cnvkit', 'delly', 'cnvpytor', 'smoove', 'mops']

    for i in range(7,13):
        truthFile = '/home/jhsun/data3/project/CNVPipe/simulation-CNVSimulator/simuGenome/sample' + str(i) + '_cnvList.bed'
        evaluateHelper(truthFile=truthFile, tools=tools, fold='10x', outputFile=out, sampleID=str(i))
    # for i in range(1,37):
    #     truthFile = '/home/jhsun/data3/project/CNVPipe/simulation-CNVSimulator/simuGenome/sample' + str(i) + '_cnvList.bed'
    #     if 1 <= i <= 6:
    #         evaluateHelper(truthFile=truthFile, tools=tools, fold='1x', outputFile=out, sampleID=str(i))
    #     elif 7 <= i <= 12:
    #         evaluateHelper(truthFile=truthFile, tools=tools, fold='10x', outputFile=out, sampleID=str(i))
    #     elif 13 <= i <= 18:
    #         evaluateHelper(truthFile=truthFile, tools=tools, fold='30x', outputFile=out, sampleID=str(i))
    #     elif 19 <= i <= 24:
    #         evaluateHelper(truthFile=truthFile, tools=tools, fold='0.1x', outputFile=out, sampleID=str(i))
    #     elif 25 <= i <= 30:
    #         evaluateHelper(truthFile=truthFile, tools=tools, fold='0.5x', outputFile=out, sampleID=str(i))
    #     else:
    #         evaluateHelper(truthFile=truthFile, tools=tools, fold='5x', outputFile=out, sampleID=str(i))

