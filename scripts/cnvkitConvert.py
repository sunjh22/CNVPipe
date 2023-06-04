#! /usr/bin/env python

# Merge consecutive bins in CNVKit results

from utils import readCNVFile, mergeConsecutiveSegments

inputFile = snakemake.input[0]
outputFile = snakemake.output[0]

cnvList = mergeConsecutiveSegments(readCNVFile(inputFile, tool='CNVKit'))

with open(outputFile, 'w') as f:
    for x in cnvList:
        print(*x, sep='\t', file=f)
