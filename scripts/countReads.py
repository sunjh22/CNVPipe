# /usr/bin/env python

import sys
import pysam
import statistics
import os

if len(sys.argv) != 5:
    sys.exit("Wrong input files! Usage: python countReads.py bamFile boundaryFile countsFile medianCountsFile")

inputBam = sys.argv[1]
inputBins = sys.argv[2]
outputCount = sys.argv[3]
outputMedianStat = sys.argv[4]

samplename = os.path.basename(inputBam).split('.')[0]

bamfile = pysam.AlignmentFile(inputBam, 'rb')

with open(inputBins, 'r') as f:
    bins = [x.strip().split('\t')[:4] for x in f if not x.startswith('chrom')]

count = []
for bin in bins:
    tmp_count = 0
    for read in bamfile.fetch(bin[0], start=int(bin[1]), end=int(bin[3])):
        tmp_count += 1
    count.append(tmp_count)

countSum = sum(count)
countMedian = statistics.median(count)

with open(outputCount, 'w') as f:
    print('chrom\tstart\tabsStart\tcount', file=f)
    for i in range(len(bins)):
        print('\t'.join(bins[i][:3]), count[i], sep='\t', file=f)

with open(outputMedianStat, 'a') as f:
    print(samplename, countSum, countMedian, sep='\t', file=f)