#!/usr/bin/env python3
"""Convert Delly CNV calls into standard BED format, resolving conflicting CNVs."""

from utils import readCNVFile, resolveConflictCNVs

inputFile = snakemake.input[0]
outputFile = snakemake.output[0]

cnvList = resolveConflictCNVs(readCNVFile(inputFile, tool='Delly'))

with open(outputFile, 'w') as f:
    for x in cnvList:
        print(*x, sep='\t', file=f)
