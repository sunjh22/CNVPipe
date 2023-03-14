#! /usr/bin/env python
# Filter confilct or complex regions for Lumpy results

import sys

def overlap(cnv1, cnv2):
    '''Calculate overlap proportion of two CNVs, return True if it is larger than some threshold
    - cnv1: ground truth CNV
    - cnv2: detected CNV
    - return: size of overlapped region between two CNVs
    '''

    c1, s1, e1 = cnv1[0], int(cnv1[1]), int(cnv1[2])
    c2, s2, e2 = cnv2[0], int(cnv2[1]), int(cnv2[2])
    cnvLen1, cnvLen2 = e1 - s1, e2 - s2
    overlap = 0
    assert cnvLen1 > 0 and cnvLen2 > 0, "The length of CNV is < 0, something wrong!"
    if c1 == c2:
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

    return overlap

if __name__ == '__main__':

    inputFile = sys.argv[1]     # smoove temp CNV bed file, say temp/smoove/sample1.bed
    outputFile = sys.argv[2]

    cnvs = []
    with open(inputFile, 'r') as f:
        for line in f:
            if line.find('_') != -1:
                continue
            line = line.strip().split('\t')
            cnvs.append(line)

    out = open(outputFile, 'w')
    for x in cnvs:
        cnvLen = int(x[2]) - int(x[1])
        flag = 0
        for y in cnvs:
            if x !=y:
                if overlap(x[:3], y[:3]) > 0:
                    flag = 1

        # the length of identified CNV should larger than 1kb
        if flag == 0 and cnvLen >= 1000:
            print(*x, sep='\t', file=out)

    out.close()