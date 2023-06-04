#! /usr/bin/env python

# Filter confilct or complex regions for Lumpy results

from utils import readCNVFile, resolveConflictCNVs

inputFile = snakemake.input[0]
outputFile = snakemake.output[0]

cnvList = resolveConflictCNVs(readCNVFile(inputFile, tool='Delly'))

with open(outputFile, 'w') as f:
    for x in cnvList:
        print(*x, sep='\t', file=f)
