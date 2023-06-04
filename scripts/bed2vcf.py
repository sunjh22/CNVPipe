#! /usr/bin/env python
# Transform merged CNV list from bed format into vcf format for duphold analysis

import sys
from datetime import date
import os

def readChrLength(infile):
    chrLength = {}
    with open(infile, 'r') as f:
        for x in f:
            if x.find('_') != -1 and not x.startswith('NC'):
                continue
            x = x.strip().split('\t')
            chrLength[x[0]] = int(x[1])
    
    return chrLength


if __name__ == "__main__":
    
    inputBedFile = sys.argv[1]      # a CNV bed file after merging, say test.bed
    inputFaiFile = sys.argv[2]    # fasta index, say hg38.fa.fai
    outputVcfFile = sys.argv[3]     # output CNV vcf file, say test.vcf

    Out = open(outputVcfFile, 'w')

    today = date.today().strftime("%Y%m%d")
    chrLength = readChrLength(inputFaiFile)
    sample = os.path.basename(inputBedFile).split('.')[0]

    print('##fileformat=VCF4.2', file=Out)
    print('##fileDate={:s}'.format(today), file=Out)
    print('##FILTER=<ID=PASS,Description="All filters passed">', file=Out)
    print('##ALT=<ID=DUP,Description="Copy number gain">', file=Out)
    print('##ALT=<ID=DEL,Description="Copy number loss">', file=Out)
    print('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of CNV">', file=Out)
    print('##INFO=<ID=TN,Number=1,Type=Integer,Description="Number of tools overlapped \
        for this CNV">', file=Out)
    print('##INFO=<ID=TNa,Number=1,Type=String,Description="Name of tools overlapped \
        for this CNV">', file=Out)
    print('##INFO=<ID=SAMPLE,Number=1,Type=String,Description="Sample name">', file=Out)
    print('##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Integer copy number">', file=Out)
    print('##FORMAT=<ID=AS,Number=1,Type=Float,Description="Accumulated overlapped \
        fraction/score, the bigger this value, the more confidence of CNV">', file=Out)
    print('##reference=/data/jinwf/jhsun/refs/hg38/analysisSet/hg38.analysisSet.fa', file=Out)

    for key in chrLength.keys():
        print('##contig=<ID={:s},length={:d}>'.format(key, chrLength[key]), file=Out)

    print('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT', sample, sep='\t', file=Out)

    id = 0
    with open(inputBedFile, 'r') as f:
        for x in f:
            if x.startswith('chromosome'):
                continue
            x = x.strip().split('\t')
            chrom = x[0]
            pos = int(x[1])
            end = int(x[2])     # bed to vcf, end should minus 1
            cn = int(x[3])
            tools = x[4]
            toolNum = int(x[5])
            accumScore = float(x[6])
            spl = x[7]
            id += 1
            ref = 'N'
            alt = 'DUP' if cn > 2 else 'DEL'
            qual = 1000
            filt = 'PASS'
            print(chrom, pos, '{:0>3d}'.format(id), ref, alt, qual, filt, 
            'END={:d};TNa={:s};TN={:d};SAMPLE={:s}'.format(end-1, tools, toolNum, spl), 'CN:AS', 
            "{:d}:{:.1f}".format(cn, accumScore), sep='\t', file=Out)

