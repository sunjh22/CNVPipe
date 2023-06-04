#! /usr/bin/env python

# Merge CNV calling results from 5 tools, assign 'accumulative Score' (1. AS).
# Only for high-depth data, we use very strict overlap threshold, that is, reciprocal overlap over
# 0.75. For median-depth and low-depth data, we use a looser one: the max overlap between two CNVs
# over 0.3 will be counted as a true overlap, thus breakpoints will be extended.

import os
from utils import readCNVFile, mergeCNVFromTools

sample = os.path.basename(snakemake.input[0]).split('.')[0]
cnvs = []

if snakemake.params[0] == 'high':
    smoove = snakemake.input[0]
    delly = snakemake.input[1]
    cnvkit = snakemake.input[2]
    cnvpytor = snakemake.input[3]
    mops = snakemake.input[4]
    outputFile = snakemake.output[0]
    # there is priority for keeping CNVs when merging
    cnvfiles = [smoove, delly, cnvkit, cnvpytor, mops]
    print('Merging CNV results from: ', cnvfiles)
    cnvtools = ['smoove', 'delly', 'cnvkit', 'cnvpytor', 'mops']
    for i, cnvfile in enumerate(cnvfiles):
        for cnv in readCNVFile(cnvfile, 'merge'):
            cnv.append(cnvtools[i])
            # five columns in cnvs: chr, start, end, cn, tool
            cnvs.append(cnv)
elif snakemake.params[0] == 'median':
    mops = snakemake.input[0]
    cnvkit = snakemake.input[1]
    delly = snakemake.input[2]
    cnvpytor = snakemake.input[3]
    smoove = snakemake.input[4]
    outputFile = snakemake.output[0]
    cnvfiles = [mops, cnvkit, delly, cnvpytor, smoove]
    print('Merging CNV results from: ', cnvfiles)
    cnvtools = ['mops', 'cnvkit', 'delly', 'cnvpytor', 'smoove']
    for i, cnvfile in enumerate(cnvfiles):
        for cnv in readCNVFile(cnvfile, 'merge'):
            cnv.append(cnvtools[i])
            cnvs.append(cnv)
elif snakemake.params[0] == 'low':
    mops = snakemake.input[0]
    cnvkit = snakemake.input[1]
    cnvpytor = snakemake.input[2]
    outputFile = snakemake.output[0]
    cnvfiles = [mops, cnvkit, cnvpytor]
    print('Merging CNV results from: ', cnvfiles)
    cnvtools = ['mops', 'cnvkit', 'cnvpytor']
    for i, cnvfile in enumerate(cnvfiles):
        for cnv in readCNVFile(cnvfile, 'merge'):
            cnv.append(cnvtools[i])
            cnvs.append(cnv)

mergedCnvs = mergeCNVFromTools(cnvs, min_threshold=0.75, max_threshold=0.95)

with open(outputFile, 'w') as f:
    print('chromosome', 'start', 'end', 'cn', 'tools', 'toolNum', 'accumScore', 'sample', sep='\t', file=f)
    for cnv in mergedCnvs:
        print(*cnv, sample, sep='\t', file=f)

