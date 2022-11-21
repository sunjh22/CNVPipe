#! /usr/bin/env python
# Filter confilct or complex regions for Lumpy results

import sys
from mergeCNV import overlapLen

inputFile = sys.argv[1]     # smoove temp CNV bed file, say temp/smoove/sample1.bed
outputFile = sys.argv[2]

cnvs = []
with open(inputFile, 'r') as f:
    for line in f:
        line = line.strip().split('\t')
        cnvs.append(line)

out = open(outputFile, 'w')
for x in cnvs:
    cnvLen = int(x[2]) - int(x[1])
    flag = 0
    for y in cnvs:
        if x !=y:
            if overlapLen(x[:3], y[:3]) > 0:
                flag = 1

    # the length of identified CNV should larger than 1kb
    if flag == 0 and cnvLen >= 1000:
        print(*x, sep='\t', file=out)

out.close()