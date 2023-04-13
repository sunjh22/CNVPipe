#! /urs/bin/env python
# Usage: python findFPCNV.py

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
    if min(cnvProp1, cnvProp2) > 0.3:
    # if cnvProp1 >= 0.3:
        return True
    else:
        return False


def evaluate(truthFile, callFile, outputFile1, outputFile2, dup=True):
    '''Calculate sensitivity and FDR for single tool and merged results. We generally take CNVs
    detected by more than one tool in merged results as input.
    - truthFile: a file with ground truth CNVs
    - callFile: a file with detected CNVs
    - return: sensitivity and FDR
    '''

    # define some threshold to filter CNV for merged result
    dupholdScoreThe = 90    # cannot be adjusted
    toolNumThe = 4         # CNVs called by at least how many tools, equal to accumScoreThe > 0

    # read detected CNVs, read more columns for CNV filtering for merged CNV set
    callCnvs = []
    with open(callFile, 'r') as f:
        for line in f:
            if line.startswith('chrom'):
                continue
            x = line.strip().split('\t')
            dupholdScore = float(x[6])
            toolName = x[7].split(',')
            toolNum = int(x[8])
            cnvfilter = x[9]
            cn = int(x[3])
            if not dup and cn > 2:
                continue
            callCnvs.append(x)
            # if (toolNum >= toolNumThe or 'smoove' in toolName or 'delly' in toolName) and dupholdScore > 0:
            # if (toolNum >= toolNumThe or 'smoove' in toolName or 'delly' in toolName) and cnvfilter == 'True':
                # callCnvs.append(x)

    # read ground truth CNVs
    truthCnvs = []
    with open(truthFile, 'r') as f:
        for line in f:
            if line.startswith('chrom'):
                continue
            x = line.strip().split('\t')[:4]
            if not dup:
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

    # the reason for tp and observTP is that several called CNVs might overlap with the same truth CNV.
    tp = []     # CNVs in truth set that can be identified by calling
    observTP = []       # CNVs in called set that have overlap with truth-set CNVs
    callCnvs2 = callCnvs[:]
    callCnvs3 = callCnvs[:]
    while truthCnvs:
        cnv1 = truthCnvs.pop(0)
        count = 0       # count the number of CNV in calling set that has overlap with truth CNV
        flag = 0        # mark whether the CNV in truth set has been called
        for i, x in enumerate(callCnvs):
            cnv2 = x[:4]
            if overlap(cnv1, cnv2):
                observTP.append(callCnvs2.pop(i-count))
                count += 1
                flag = 1

        if flag == 1:
            tp.append(cnv1)
        callCnvs = callCnvs2[:]

    with open(outputFile1, 'w') as f, open(outputFile2, 'w') as g:
        for x in callCnvs3:
            if x in observTP:
                print(*x, sep='\t', file=f)
            else:
                print(*x, sep='\t', file=g)


if __name__ == "__main__":

    samples = ['NA12878-1', 'NA12878-2', 'CHM13', 'AK1', 'HG002', 'HG00514', 'HG00733', 'NA19240']

    for sampleID in samples:
        outputFile1 = '/home/jhsun/data3/project/CNVPipe/realAnalysis-10x/evaluation/' + sampleID + '.observTP.bed'
        outputFile2 = '/home/jhsun/data3/project/CNVPipe/realAnalysis-10x/evaluation/' + sampleID + '.FP.bed'
        truthFile = '/home/jhsun/data3/project/CNVPipe/realAnalysis/truthSet/' + sampleID + '-SVset.1kb.bed'
        callFile = '/home/jhsun/data3/project/CNVPipe/realAnalysis-10x/res/merge/' + sampleID + '.bed'
        if sampleID in ['sample13', 'sample14']:
            evaluate(truthFile=truthFile, callFile=callFile, outputFile1=outputFile1, outputFile2=outputFile2, dup=True)
        else:
            evaluate(truthFile=truthFile, callFile=callFile, outputFile1=outputFile1, outputFile2=outputFile2, dup=False)
