#! /usr/bin/env python
# Find recurrent CNVs among all samples

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
    # cnvProp = min(cnvProp1, cnvProp2)
    cnvProp = max(cnvProp1, cnvProp2)
    return overlap, cnvProp


def getCNVs(workDir):
    new_cnvs = []
    files = os.listdir(workDir + 'res/CNVPipe/')
    for tmpFile in files:
        with open('res/CNVPipe/' + tmpFile, 'r') as f:
            sample = os.path.basename(tmpFile).split('.')[0]
            for cnv in f:
                if cnv.startswith('chrom'):
                        continue
                cnv = cnv.strip().split('\t')
                chrom = cnv[0]
                start = int(cnv[1])
                end = int(cnv[2])
                cn = int(cnv[3])
                dupScore = int(cnv[6])
                toolNum = int(cnv[8])
                cnvfilter = cnv[9]
                goodScore = int(cnv[10])
                if toolNum >=3 and cnvfilter == 'True' and dupScore == 100 and goodScore == 100:
                    new_cnv = [chrom, start, end, cn, sample]
                    new_cnvs.append(new_cnv)

    return new_cnvs


def main(outputFile, workDir, controlNames, sampleNumThe):
    cnvs = getCNVs(workDir)

    overlapPropThreshold = 0.8
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
            if overlapProp < overlapPropThreshold:
                continue
            else:
                assert cnv1[-1] != cnv2[-1], f"Overlapped CNV from same sample {cnv1[-1]:s}! Please make sure these conflicts are solved before merging"
                tmpCnv[1:3] = [max(cnv1[1], cnv2[1]), min(cnv1[2], cnv2[2])]
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
        print('chromosome', 'start', 'end', 'cn', 'samples', 'sampleNum', 'accumScore', sep='\t', file=f)
        for cnv in mergedCnvs:
            sampleName = cnv[4]
            sampleNum = cnv[5]
            if sampleNum < int(sampleNumThe):
                continue
            for tmpSample in sampleName:
                if tmpSample in controlNames:
                    continue
            print(*cnv, sep='\t', file=f)


if __name__ == '__main__':

    outputFile = sys.argv[1]
    workDir = sys.argv[2]
    sampleNumThe = sys.argv[3]  # threshold of the number of samples to keep the recurrent CNVs
    controlNames = sys.argv[4:]  # a list of control sample names
    
    # print(controlNames)

    main(outputFile, workDir, controlNames, sampleNumThe)