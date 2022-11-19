#! /usr/bin/env python
# Merge simulated ground truth CNV and SV into one file, do some filtering

import sys

def mergeTruthCNV(cnv1, cnv2):
    c1, s1, e1, cn1 = cnv1[0], int(cnv1[1]), int(cnv1[2]), int(cnv1[3])
    c2, s2, e2, cn2 = cnv2[0], int(cnv2[1]), int(cnv2[2]), int(cnv2[3])
    if c1 == c2:
        if (cn1>2 and cn2>2) or (cn1<2 and cn2<2):
            if s1 <= s2 < e2 <= e1:
                return([c1, s1, e1, cn1])
            elif s1 <= s2 <= e1 <= e2:
                return([c1, s1, e2, cn1])
            elif s2 <= s1 < e1 <= e2:
                return([c1, s2, e2, cn1])
            elif s2 <= s1 <= e2 <= e1:
                return([c1, s2, e1, cn1])
            else:
                return(cnv1)
    
    return(cnv1)


cnvFile = sys.argv[1]   # CNV truth set
svFile = sys.argv[2]    # SV truth set
outputFile = sys.argv[3]

cnvs = []
with open(cnvFile, 'r') as f:
    for line in f:
        if line.startswith('chr'):
            continue
        x = line.strip().split('\t')
        x[0] = 'chr' + x[0]
        size = int(x[2]) - int(x[1])
        if size < 5000:
            continue
        cn = int(x[-1])
        if cn == 0:
            continue
        elif cn < 0:
            x[-1] = '1'
        else:
            x[-1] = '3'
        cnvs.append(x)
# print(cnvs)

svs = []
with open(svFile, 'r') as f:
    for line in f:
        if line.startswith('chr'):
            continue
        x = line.strip().split('\t')
        chrom = 'chr' + x[0]
        start = x[1]
        end = x[3]
        size = int(start) - int(end)
        if size < 5000:
            continue
        Type = x[-1]
        if Type == 'Deletion':
            cn = '1'
        elif Type == 'TandemDup':
            cn = '3'
        else:
            continue
        svs.append([chrom, start, end, cn])
# print(svs)

merged = []
svs2 = svs[:]
while cnvs:
    cnv1 = cnvs.pop(0)
    count = 0
    tmpCnv = []
    for i, cnv2 in enumerate(svs):
        tmpCnv = mergeTruthCNV(cnv1, cnv2)
        if tmpCnv != cnv1:
            svs2.pop(i-count)
            count += 1
            cnv1 = tmpCnv

    merged.append(cnv1)
    svs = svs2[:]

merged.extend(svs)
with open(outputFile, 'w') as f:
    for x in merged:
        print(*x, sep='\t', file=f)
